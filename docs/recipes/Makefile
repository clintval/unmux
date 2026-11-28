# Run the pixi tasks of the cookbook workflow; CACHE=dir runs it from an empty cache in dir.

PIXI ?= pixi
CACHE ?=
CONFIG = $(if $(CACHE),--config cache=$(CACHE) results=$(CACHE)/results)

.PHONY: page run check clean

page:
	$(PIXI) run page $(CONFIG)

run:
	$(PIXI) run run $(CONFIG)

check:
	$(PIXI) run check $(CONFIG)

clean:
	$(PIXI) run clean
