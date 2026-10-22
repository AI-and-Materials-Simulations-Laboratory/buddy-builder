# Build/install hooks for the nanoHUB tool pipeline (spec H5).
# Pure Python: nothing to compile. `make install` copies the package to
# ../lib/python, where app.ipynb looks first.

INSTALL_DIR ?= ../lib/python
PACKAGE = nanobuilder

all:
	@echo "Pure Python package; nothing to build."

install:
	rm -rf $(INSTALL_DIR)/$(PACKAGE)
	mkdir -p $(INSTALL_DIR)
	cp -r $(PACKAGE) $(INSTALL_DIR)/
	find $(INSTALL_DIR)/$(PACKAGE) -name "__pycache__" -type d -prune -exec rm -rf {} +

test:
	cd .. && python -m pytest tests -q

clean:
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +

distclean: clean
	rm -rf ../lib

.PHONY: all install test clean distclean
