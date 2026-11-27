//! SAM tags tied to a read's exact bases or pulses, and how each follows changes to the read.

use std::convert::identity;

use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::record_buf::data::field::Value;

/// The value type a per-base tag has; a same-named tag of another type is a different tag.
#[derive(Clone, Copy)]
pub(crate) enum Kind {
    /// Any type.
    Any,
    /// A `Z` string.
    Text,
    /// A `B` array of any subtype.
    Array,
}

/// How a per-base tag follows a read whose bases are cut, reordered or reverse-complemented.
#[derive(Clone, Copy)]
pub(crate) enum PerBase {
    /// Left out of a changed record.
    Dropped,
    /// Cut to a changed record's segments and reversed in a flipped one; aligners leave it as sequenced.
    Native,
    /// Cut to a changed record's segments and reversed in a flipped one.
    Reversed,
    /// Bases cut to a changed record's segments and reverse-complemented in a flipped one.
    Complemented,
}

/// Every per-base tag, the type it must have, and how it follows its read.
pub(crate) const PER_BASE_TAGS: [([u8; 2], Kind, PerBase); 54] = [
    (*b"MM", Kind::Any, PerBase::Dropped),
    (*b"ML", Kind::Any, PerBase::Dropped),
    (*b"MN", Kind::Any, PerBase::Dropped),
    (*b"OQ", Kind::Text, PerBase::Reversed),
    (*b"BQ", Kind::Any, PerBase::Dropped),
    (*b"E2", Kind::Text, PerBase::Complemented),
    (*b"U2", Kind::Text, PerBase::Reversed),
    (*b"CQ", Kind::Any, PerBase::Dropped),
    (*b"CS", Kind::Any, PerBase::Dropped),
    (*b"FZ", Kind::Any, PerBase::Dropped),
    (*b"PT", Kind::Any, PerBase::Dropped),
    (*b"BI", Kind::Any, PerBase::Dropped),
    (*b"BD", Kind::Any, PerBase::Dropped),
    (*b"ip", Kind::Array, PerBase::Native),
    (*b"pw", Kind::Array, PerBase::Native),
    (*b"fi", Kind::Any, PerBase::Dropped),
    (*b"fp", Kind::Any, PerBase::Dropped),
    (*b"ri", Kind::Any, PerBase::Dropped),
    (*b"rp", Kind::Any, PerBase::Dropped),
    (*b"dq", Kind::Any, PerBase::Dropped),
    (*b"dt", Kind::Any, PerBase::Dropped),
    (*b"iq", Kind::Any, PerBase::Dropped),
    (*b"mq", Kind::Any, PerBase::Dropped),
    (*b"sq", Kind::Any, PerBase::Dropped),
    (*b"sa", Kind::Any, PerBase::Dropped),
    (*b"sx", Kind::Any, PerBase::Dropped),
    (*b"pc", Kind::Any, PerBase::Dropped),
    (*b"pq", Kind::Any, PerBase::Dropped),
    (*b"pg", Kind::Any, PerBase::Dropped),
    (*b"pa", Kind::Any, PerBase::Dropped),
    (*b"pm", Kind::Any, PerBase::Dropped),
    (*b"pd", Kind::Any, PerBase::Dropped),
    (*b"px", Kind::Any, PerBase::Dropped),
    (*b"pe", Kind::Any, PerBase::Dropped),
    (*b"sf", Kind::Any, PerBase::Dropped),
    (*b"bs", Kind::Any, PerBase::Dropped),
    (*b"pb", Kind::Any, PerBase::Dropped),
    (*b"mv", Kind::Any, PerBase::Dropped),
    (*b"cd", Kind::Array, PerBase::Reversed),
    (*b"ce", Kind::Array, PerBase::Reversed),
    (*b"ad", Kind::Array, PerBase::Reversed),
    (*b"ae", Kind::Array, PerBase::Reversed),
    (*b"bd", Kind::Array, PerBase::Reversed),
    (*b"be", Kind::Array, PerBase::Reversed),
    (*b"aq", Kind::Text, PerBase::Reversed),
    (*b"bq", Kind::Text, PerBase::Reversed),
    (*b"ac", Kind::Text, PerBase::Complemented),
    (*b"bc", Kind::Text, PerBase::Complemented),
    (*b"cu", Kind::Array, PerBase::Reversed),
    (*b"ct", Kind::Array, PerBase::Reversed),
    (*b"au", Kind::Array, PerBase::Reversed),
    (*b"at", Kind::Array, PerBase::Reversed),
    (*b"bu", Kind::Array, PerBase::Reversed),
    (*b"bt", Kind::Array, PerBase::Reversed),
];

const _: () = assert!(PER_BASE_TAGS.len() <= u64::BITS as usize);

/// The [`PER_BASE_TAGS`] index and handling of a tag whose value has its per-base type.
pub(crate) fn per_base(tag: Tag, value: &Value) -> Option<(usize, PerBase)> {
    let index = PER_BASE_TAGS
        .iter()
        .position(|(name, _, _)| tag.as_ref() == name)?;
    let (_, kind, per_base) = PER_BASE_TAGS[index];
    let typed = match kind {
        Kind::Any => true,
        Kind::Text => matches!(value, Value::String(_)),
        Kind::Array => matches!(value, Value::Array(_)),
    };
    typed.then_some((index, per_base))
}

/// A record segment: bases `[start, end)` of its read, reverse-complemented when `flipped`.
#[derive(Clone, Copy)]
pub(crate) struct Slice {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) flipped: bool,
}

/// The values of a `len`-base read cut to `slices`, reversing each flipped one through `flip`.
fn cut<T: Copy>(
    values: &[T],
    slices: &[Slice],
    len: usize,
    flip: impl Fn(T) -> T,
) -> Option<Vec<T>> {
    if values.len() != len {
        return None;
    }
    let mut cut = Vec::with_capacity(slices.iter().map(|slice| slice.end - slice.start).sum());
    for slice in slices {
        let part = values.get(slice.start..slice.end)?;
        match slice.flipped {
            true => cut.extend(part.iter().rev().map(|&value| flip(value))),
            false => cut.extend_from_slice(part),
        }
    }
    Some(cut)
}

/// A per-base `value` of a `len`-base read cut to a record's `slices`, if it has `len` entries.
pub(crate) fn follow(
    value: &Value,
    per_base: PerBase,
    slices: &[Slice],
    len: usize,
) -> Option<Value> {
    Some(match value {
        Value::String(text) => {
            let flip: fn(u8) -> u8 = match per_base {
                PerBase::Complemented => crate::iupac::complement,
                _ => identity,
            };
            Value::String(cut(&text[..], slices, len, flip)?.into())
        }
        Value::Array(array) => Value::Array(match array {
            Array::Int8(values) => Array::Int8(cut(values, slices, len, identity)?),
            Array::UInt8(values) => Array::UInt8(cut(values, slices, len, identity)?),
            Array::Int16(values) => Array::Int16(cut(values, slices, len, identity)?),
            Array::UInt16(values) => Array::UInt16(cut(values, slices, len, identity)?),
            Array::Int32(values) => Array::Int32(cut(values, slices, len, identity)?),
            Array::UInt32(values) => Array::UInt32(cut(values, slices, len, identity)?),
            Array::Float(values) => Array::Float(cut(values, slices, len, identity)?),
        }),
        _ => return None,
    })
}

/// The number of entries of a text or array value.
fn entries(value: &Value) -> Option<usize> {
    Some(match value {
        Value::String(text) => text.len(),
        Value::Array(Array::Int8(values)) => values.len(),
        Value::Array(Array::UInt8(values)) => values.len(),
        Value::Array(Array::Int16(values)) => values.len(),
        Value::Array(Array::UInt16(values)) => values.len(),
        Value::Array(Array::Int32(values)) => values.len(),
        Value::Array(Array::UInt32(values)) => values.len(),
        Value::Array(Array::Float(values)) => values.len(),
        _ => return None,
    })
}

/// A tag's value in sequencing order, for a record whose `SEQ` an aligner reverse-complemented.
pub(crate) fn sequenced(tag: Tag, value: &Value) -> Value {
    let turned = per_base(tag, value).and_then(|(_, per_base)| match per_base {
        PerBase::Reversed | PerBase::Complemented => {
            let len = entries(value)?;
            let whole = Slice {
                start: 0,
                end: len,
                flipped: true,
            };
            follow(value, per_base, &[whole], len)
        }
        PerBase::Dropped | PerBase::Native => None,
    });
    turned.unwrap_or_else(|| value.clone())
}
