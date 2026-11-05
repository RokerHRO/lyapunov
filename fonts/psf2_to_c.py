#!/usr/bin/env python3

import sys
import struct

FIRST_CHAR = 0x20
LAST_CHAR = 0x7F

PSF2_MAGIC = 0x864AB572
PSF2_HEADER_SIZE = 32


def extract_psf2(input_file, output_file):
    with open(input_file, "rb") as f:
        data = f.read()

    if len(data) < PSF2_HEADER_SIZE:
        raise ValueError("File too small for PSF2 header")

    (
        magic,
        version,
        headersize,
        flags,
        numglyph,
        bytesperglyph,
        height,
        width
    ) = struct.unpack("<8I", data[:32])

    if magic != PSF2_MAGIC:
        raise ValueError("Invalid PSF2 format")

    if width < 1 or width > 16:
        raise ValueError(
            f"uint16_t supports a maximum width of 16 pixels, "
            f"found {width}"
        )

    if height < 1:
        raise ValueError("Invalid glyph height")

    if numglyph <= LAST_CHAR:
        raise ValueError("Font does not contain enough glyphs")

    bytes_per_row = (width + 7) // 8
    required_bytes = bytes_per_row * height

    if bytesperglyph < required_bytes:
        raise ValueError("Glyph size in header is too small")

    if headersize + numglyph * bytesperglyph > len(data):
        raise ValueError("Incomplete glyph data")

    glyphs = []

    for char in range(FIRST_CHAR, LAST_CHAR + 1):
        offset = headersize + char * bytesperglyph
        glyph_data = data[offset:offset + bytesperglyph]

        rows = []

        for y in range(height):
            row_offset = y * bytes_per_row
            row_bytes = glyph_data[
                row_offset:row_offset + bytes_per_row
            ]

            # Convert MSB-first pixel data to uint16_t
            value = int.from_bytes(row_bytes, "big")
            value >>= bytes_per_row * 8 - width

            rows.append(value)

        glyphs.append(rows)

    with open(output_file, "w", encoding="ascii") as f:
        f.write("#include <stdint.h>\n\n")
        f.write(f"#define FONT_FIRST_CHAR 0x{FIRST_CHAR:02X}\n")
        f.write(f"#define FONT_LAST_CHAR  0x{LAST_CHAR:02X}\n")
        f.write(f"#define FONT_WIDTH      {width}\n")
        f.write(f"#define FONT_HEIGHT     {height}\n\n")

        f.write(f"const uint16_t font[96][{height}] = {{\n")

        for i, rows in enumerate(glyphs):
            char = FIRST_CHAR + i

            f.write("    {\n")

            for start in range(0, len(rows), 8):
                values = rows[start:start + 8]
                f.write("        ")
                f.write(", ".join(f"0x{v:04X}" for v in values))

                if start + 8 < len(rows):
                    f.write(",")

                f.write("\n")

            f.write(f"    }},  /* 0x{char:02X} */\n")

        f.write("};\n")

    print(f"Generated {output_file} successfully.")
    print(f"Glyphs: {len(glyphs)}")
    print(f"Size:   {width} x {height}")
    print(f"Bytes:  {bytesperglyph} per glyph")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.psf output.c")
        sys.exit(1)

    extract_psf2(sys.argv[1], sys.argv[2])

