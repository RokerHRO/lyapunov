#!/usr/bin/python3

import urllib.request
import sys
from pathlib import Path

URL = sys.argv[1]

data = Path(URL).read_bytes()

print("Read")

# PSF1-Header ueberpruefen
assert data[0:2] == bytes([0x36, 0x04]), "Ungueltiges PSF1-Format"
assert data[3] == 16, "Unerwartete Glyphengroesse"

glyph_size = data[3]
first_char = 0x20
last_char = 0x7F
count = last_char - first_char + 1

# Glyphendaten extrahieren
offset = 4 + first_char * glyph_size
glyphs = data[offset:offset + count * glyph_size]

assert len(glyphs) == 1536, "Unvollstaendige Glyphendaten"

# C-Quelldatei generieren
with open(URL + ".c", "w") as f:
    f.write('#include <stdint.h>\n\n')
    f.write('#define FONT_FIRST_CHAR 0x20\n')
    f.write('#define FONT_LAST_CHAR  0x7F\n')
    f.write('#define FONT_WIDTH      8\n')
    f.write('#define FONT_HEIGHT     16\n\n')

    f.write('const uint8_t zap_light16[96][16] = {\n')

    for i in range(count):
        glyph = glyphs[i * 16:(i + 1) * 16]

        f.write('    { ')
        f.write(', '.join(f'0x{b:02X}' for b in glyph))
        f.write(' },  /* 0x%02X */\n' % (first_char + i))

    f.write('};\n')

print(URL + ".c erfolgreich erzeugt.")
