RIS is a simple method for detecting and validating chimeric reads from mappings to multiple genomes.
=========================================================================

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://opensource.org/licenses/GPL-3.0
    :alt: GPLv3 License

.. contents::
   :local:
   :depth: 2

Introduction
^^^^^^^^^^^^
RIS is a simple method for detecting and validating chimeric reads from mappings to multiple genomes.

Installation
^^^^^^^^^^^^

Building From Source
""""""""""""""""""""
::

    $ git clone https://github.com/alevar/ris.git
    $ cd ris

**Requirements**
Install the Rust compiler with `rustup`

1. Install https://rustup.rs/.
2. Test installation

::

    $ rustup override set stable
    $ rustup update stable

**Building**

::

    $ cargo build --release

Usage
^^^^^

::

    $ ris -h
