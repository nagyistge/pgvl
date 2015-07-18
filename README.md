# Pretty Good Vision Library {#mainpage}

Like the name says, I hope this ends up being a pretty good vision library.
I aim to make it efficient, succinct, and well-documented.

## Prerequisites

- cmake
- libgtest-dev
- libeigen3-dev
- libsdl2-dev

## Building

Always do an out-of-source build.

    $ tar -xJf pgvl.tar.xz
    $ mkdir pgvl-build
    $ cd pgvl-build
    $ cmake ../pgvl
    $ make

To make the documentation:

    $ make doc
    $ xdg-open doc/html/index.html

To run the tests:

    $ make test
