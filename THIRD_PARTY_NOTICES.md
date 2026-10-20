# Philox / Random123

`src/colorplane_rng.f90` is an independent Fortran implementation of
Philox4x32-10, licensed with Colorplane under GPL-3.0-or-later.
Its algorithm constants and the known-answer values used in the tests are
documented by D. E. Shaw Research's [Random123 project](https://github.com/DEShawResearch/random123),
particularly [philox.h](https://github.com/DEShawResearch/random123/blob/main/include/Random123/philox.h)
and [kat_vectors](https://github.com/DEShawResearch/random123/blob/main/tests/kat_vectors).
No Random123 headers or binaries are bundled. The upstream notice follows for attribution.

Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice,
  this list of conditions, and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions, and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the name of D. E. Shaw Research nor the names of its contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
