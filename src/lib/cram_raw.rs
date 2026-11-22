//! The raw bytes of CRAM 3 containers, for the slice header fields noodles keeps private.

use std::io::{self, Read};

/// A reader that keeps a copy of every byte read since the last [`Recorder::clear`].
pub(crate) struct Recorder<R> {
    inner: R,
    bytes: Vec<u8>,
}

impl<R> Recorder<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            bytes: Vec::new(),
        }
    }

    pub(crate) fn clear(&mut self) {
        self.bytes.clear();
    }

    pub(crate) fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

impl<R: Read> Read for Recorder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes.extend_from_slice(&buf[..n]);
        Ok(n)
    }
}

const RAW: u8 = 0;
const EXTERNAL_DATA: u8 = 4;

/// The reference bases a slice embeds: the external block that holds them and the first position.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct EmbeddedReference {
    pub(crate) content_id: i32,
    pub(crate) start: usize,
}

/// One block of a container body, borrowed.
struct Block<'a> {
    method: u8,
    content_type: u8,
    content_id: i32,
    data: &'a [u8],
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

fn byte(src: &mut &[u8]) -> io::Result<u8> {
    let (&first, rest) = src
        .split_first()
        .ok_or_else(|| invalid("a CRAM block ends early"))?;
    *src = rest;
    Ok(first)
}

fn bytes<'a>(src: &mut &'a [u8], n: usize) -> io::Result<&'a [u8]> {
    if src.len() < n {
        return Err(invalid("a CRAM block ends early"));
    }
    let (head, rest) = src.split_at(n);
    *src = rest;
    Ok(head)
}

fn read_itf8(src: &mut &[u8]) -> io::Result<i32> {
    let first = byte(src)?;
    let extra = first.leading_ones().min(4) as usize;
    let mut n = u32::from(first) & (0x7f >> extra.min(3));
    for i in 0..extra {
        let b = u32::from(byte(src)?);
        n = if extra == 4 && i == 3 {
            (n << 4) | (b & 0x0f)
        } else {
            (n << 8) | b
        };
    }
    Ok(n as i32)
}

fn skip_ltf8(src: &mut &[u8]) -> io::Result<()> {
    let extra = byte(src)?.leading_ones() as usize;
    bytes(src, extra).map(|_| ())
}

pub(crate) fn write_itf8(dst: &mut Vec<u8>, n: i32) {
    let u = n as u32;
    match u {
        0..0x80 => dst.push(u as u8),
        0x80..0x4000 => dst.extend([0x80 | (u >> 8) as u8, u as u8]),
        0x4000..0x20_0000 => dst.extend([0xc0 | (u >> 16) as u8, (u >> 8) as u8, u as u8]),
        0x20_0000..0x1000_0000 => dst.extend([
            0xe0 | (u >> 24) as u8,
            (u >> 16) as u8,
            (u >> 8) as u8,
            u as u8,
        ]),
        _ => dst.extend([
            0xf0 | (u >> 28) as u8,
            (u >> 20) as u8,
            (u >> 12) as u8,
            (u >> 4) as u8,
            (u & 0x0f) as u8,
        ]),
    }
}

fn read_block<'a>(src: &mut &'a [u8]) -> io::Result<Block<'a>> {
    let method = byte(src)?;
    let content_type = byte(src)?;
    let content_id = read_itf8(src)?;
    let size = usize::try_from(read_itf8(src)?).map_err(|_| invalid("a negative block size"))?;
    read_itf8(src)?;
    let data = bytes(src, size)?;
    bytes(src, 4)?;
    Ok(Block {
        method,
        content_type,
        content_id,
        data,
    })
}

fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = flate2::Crc::new();
    crc.update(bytes);
    crc.sum()
}

/// The reference a slice embeds, read from its raw bytes, if it embeds one in a raw slice header.
pub(crate) fn embedded_reference(slice: &[u8]) -> io::Result<Option<EmbeddedReference>> {
    let mut src = slice;
    let block = read_block(&mut src)?;
    if block.method != RAW {
        return Ok(None);
    }
    let mut header = block.data;
    let reference_sequence_id = read_itf8(&mut header)?;
    let start = read_itf8(&mut header)?;
    for _ in 0..2 {
        read_itf8(&mut header)?;
    }
    skip_ltf8(&mut header)?;
    read_itf8(&mut header)?;
    for _ in 0..read_itf8(&mut header)? {
        read_itf8(&mut header)?;
    }
    let content_id = read_itf8(&mut header)?;
    Ok(match (reference_sequence_id, usize::try_from(start)) {
        (0.., Ok(start)) if content_id >= 0 => Some(EmbeddedReference { content_id, start }),
        _ => None,
    })
}

/// Split container `raw` into its header and its body, whose length the header opens with.
pub(crate) fn split(raw: &[u8]) -> io::Result<(&[u8], &[u8])> {
    let length = raw
        .first_chunk::<4>()
        .and_then(|length| usize::try_from(i32::from_le_bytes(*length)).ok())
        .ok_or_else(|| invalid("a CRAM container header ends early"))?;
    let at = raw
        .len()
        .checked_sub(length)
        .ok_or_else(|| invalid("a CRAM container is shorter than its body"))?;
    Ok(raw.split_at(at))
}

/// The raw bytes of the slice at `landmark` in container `raw`.
pub(crate) fn slice<'a>(
    raw: &'a [u8],
    landmarks: &[usize],
    landmark: usize,
) -> io::Result<&'a [u8]> {
    let (_, body) = split(raw)?;
    let start = *landmarks
        .get(landmark)
        .ok_or_else(|| invalid("a CRAM slice has no landmark"))?;
    let end = landmarks.get(landmark + 1).copied().unwrap_or(body.len());
    body.get(start..end)
        .ok_or_else(|| invalid("a CRAM landmark is out of range"))
}

/// Container `raw` rebuilt so the slice at `landmark` embeds `bases` as its reference.
pub(crate) fn with_embedded_bases(
    raw: &[u8],
    landmarks: &[usize],
    landmark: usize,
    embedded: EmbeddedReference,
    bases: &[u8],
) -> io::Result<Vec<u8>> {
    let (header, body) = split(raw)?;
    let start = *landmarks
        .get(landmark)
        .ok_or_else(|| invalid("a CRAM slice has no landmark"))?;
    let end = landmarks.get(landmark + 1).copied().unwrap_or(body.len());
    let mut src = body
        .get(start..end)
        .ok_or_else(|| invalid("a CRAM landmark is out of range"))?;
    let (at, len) = loop {
        let before = src.len();
        let block = read_block(&mut src)?;
        if block.content_type == EXTERNAL_DATA && block.content_id == embedded.content_id {
            break (end - before, before - src.len());
        }
        if src.is_empty() {
            return Err(invalid("a CRAM slice lacks its embedded reference block"));
        }
    };

    let mut block = vec![RAW, EXTERNAL_DATA];
    write_itf8(&mut block, embedded.content_id);
    let size =
        i32::try_from(bases.len()).map_err(|_| invalid("an embedded reference is too long"))?;
    write_itf8(&mut block, size);
    write_itf8(&mut block, size);
    block.extend_from_slice(bases);
    block.extend(crc32(&block).to_le_bytes());

    let mut new_body = Vec::with_capacity(body.len() + block.len());
    new_body.extend_from_slice(&body[..at]);
    new_body.extend_from_slice(&block);
    new_body.extend_from_slice(&body[at + len..]);
    let grown = |offset: usize| {
        if offset > start {
            offset + block.len() - len
        } else {
            offset
        }
    };

    let mut fields = header
        .get(4..)
        .ok_or_else(|| invalid("a CRAM container header ends early"))?;
    let total = fields.len();
    for _ in 0..4 {
        read_itf8(&mut fields)?;
    }
    skip_ltf8(&mut fields)?;
    skip_ltf8(&mut fields)?;
    read_itf8(&mut fields)?;
    let kept = &header[4..4 + total - fields.len()];

    let body_len =
        i32::try_from(new_body.len()).map_err(|_| invalid("a CRAM container is too long"))?;
    let mut out = body_len.to_le_bytes().to_vec();
    out.extend_from_slice(kept);
    write_itf8(&mut out, landmarks.len() as i32);
    for &offset in landmarks {
        write_itf8(&mut out, grown(offset) as i32);
    }
    out.extend(crc32(&out).to_le_bytes());
    out.extend_from_slice(&new_body);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_itf8_round_trips() {
        for n in [
            0,
            1,
            127,
            128,
            16_383,
            16_384,
            2_097_151,
            2_097_152,
            268_435_455,
            268_435_456,
            i32::MAX,
            -1,
            -2,
        ] {
            let mut buf = Vec::new();
            write_itf8(&mut buf, n);
            let mut src = &buf[..];
            assert_eq!(read_itf8(&mut src).unwrap(), n);
            assert!(src.is_empty());
        }
    }

    #[test]
    fn test_skip_ltf8_skips_by_the_leading_ones() {
        for (bytes, len) in [(&[0x05][..], 1), (&[0x80, 0][..], 2), (&[0xff; 9][..], 9)] {
            let mut src = bytes;
            skip_ltf8(&mut src).unwrap();
            assert_eq!(bytes.len() - src.len(), len);
        }
    }
}
