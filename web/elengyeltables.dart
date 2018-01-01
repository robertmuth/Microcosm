library elengyel;

import 'dart:typed_data';

/*
 * Based on http://www.terathon.com/voxels/
 */

List<int> _Flatten(List<List<int>> lst, int rowLength) {
  List<int> out = <int>[];
  for (List<int> l in lst) {
    assert(l.length <= rowLength);
    out.addAll(l);
    for (int i = l.length; i < rowLength; ++i) out.add(0);
  }
  assert(out.length == rowLength * lst.length);
  return out;
}

final List<int> _hotCornerMapOriginal = <int>[
  // this comment prevents reformatting of this table
  0x00, 0x01, 0x01, 0x03, 0x01, 0x03, 0x02, 0x04,
  0x01, 0x02, 0x03, 0x04, 0x03, 0x04, 0x04, 0x03,
  0x01, 0x03, 0x02, 0x04, 0x02, 0x04, 0x06, 0x0C,
  0x02, 0x05, 0x05, 0x0B, 0x05, 0x0A, 0x07, 0x04,
  0x01, 0x02, 0x03, 0x04, 0x02, 0x05, 0x05, 0x0A,
  0x02, 0x06, 0x04, 0x0C, 0x05, 0x07, 0x0B, 0x04,
  0x03, 0x04, 0x04, 0x03, 0x05, 0x0B, 0x07, 0x04,
  0x05, 0x07, 0x0A, 0x04, 0x08, 0x0E, 0x0E, 0x03,
  0x01, 0x02, 0x02, 0x05, 0x03, 0x04, 0x05, 0x0B,
  0x02, 0x06, 0x05, 0x07, 0x04, 0x0C, 0x0A, 0x04,
  0x03, 0x04, 0x05, 0x0A, 0x04, 0x03, 0x07, 0x04,
  0x05, 0x07, 0x08, 0x0E, 0x0B, 0x04, 0x0E, 0x03,
  0x02, 0x06, 0x05, 0x07, 0x05, 0x07, 0x08, 0x0E,
  0x06, 0x09, 0x07, 0x0F, 0x07, 0x0F, 0x0E, 0x0D,
  0x04, 0x0C, 0x0B, 0x04, 0x0A, 0x04, 0x0E, 0x03,
  0x07, 0x0F, 0x0E, 0x0D, 0x0E, 0x0D, 0x02, 0x01,
  0x01, 0x02, 0x02, 0x05, 0x02, 0x05, 0x06, 0x07,
  0x03, 0x05, 0x04, 0x0A, 0x04, 0x0B, 0x0C, 0x04,
  0x02, 0x05, 0x06, 0x07, 0x06, 0x07, 0x09, 0x0F,
  0x05, 0x08, 0x07, 0x0E, 0x07, 0x0E, 0x0F, 0x0D,
  0x03, 0x05, 0x04, 0x0B, 0x05, 0x08, 0x07, 0x0E,
  0x04, 0x07, 0x03, 0x04, 0x0A, 0x0E, 0x04, 0x03,
  0x04, 0x0A, 0x0C, 0x04, 0x07, 0x0E, 0x0F, 0x0D,
  0x0B, 0x0E, 0x04, 0x03, 0x0E, 0x02, 0x0D, 0x01,
  0x03, 0x05, 0x05, 0x08, 0x04, 0x0A, 0x07, 0x0E,
  0x04, 0x07, 0x0B, 0x0E, 0x03, 0x04, 0x04, 0x03,
  0x04, 0x0B, 0x07, 0x0E, 0x0C, 0x04, 0x0F, 0x0D,
  0x0A, 0x0E, 0x0E, 0x02, 0x04, 0x03, 0x0D, 0x01,
  0x04, 0x07, 0x0A, 0x0E, 0x0B, 0x0E, 0x0E, 0x02,
  0x0C, 0x0F, 0x04, 0x0D, 0x04, 0x0D, 0x03, 0x01,
  0x03, 0x04, 0x04, 0x03, 0x04, 0x03, 0x0D, 0x01,
  0x04, 0x0D, 0x03, 0x01, 0x03, 0x01, 0x01, 0x00
];

final Uint8List hotCornerMap = new Uint8List.fromList(_hotCornerMapOriginal);

List<List<int>> _configClassOriginal = [
  // vertex count, triangle count, vertex indices
  [0, 0],
  [3, 1, 0, 1, 2],
  [6, 2, 0, 1, 2, 3, 4, 5],
  [4, 2, 0, 1, 2, 0, 2, 3],
  [5, 3, 0, 1, 4, 1, 3, 4, 1, 2, 3],
  [7, 3, 0, 1, 2, 0, 2, 3, 4, 5, 6],
  [9, 3, 0, 1, 2, 3, 4, 5, 6, 7, 8],
  [8, 4, 0, 1, 4, 1, 3, 4, 1, 2, 3, 5, 6, 7],
  [8, 4, 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7],
  [12, 4, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
  [6, 4, 0, 4, 5, 0, 1, 4, 1, 3, 4, 1, 2, 3],
  [6, 4, 0, 5, 4, 0, 4, 1, 1, 4, 3, 1, 3, 2],
  [6, 4, 0, 4, 5, 0, 3, 4, 0, 1, 3, 1, 2, 3],
  [6, 4, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5],
  [7, 5, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 5, 6],
  [9, 5, 0, 4, 5, 0, 3, 4, 0, 1, 3, 1, 2, 3, 6, 7, 8],
];

const int configClassRowLength = 17;

final Uint8List configClass = new Uint8List.fromList(
    _Flatten(_configClassOriginal, configClassRowLength));

final List<List<int>> _vertexDataOriginal = [
  [],
  [0x6201, 0x5102, 0x3304],
  [0x6201, 0x2315, 0x4113],
  [0x5102, 0x3304, 0x2315, 0x4113],
  [0x5102, 0x4223, 0x1326],
  [0x3304, 0x6201, 0x4223, 0x1326],
  [0x6201, 0x2315, 0x4113, 0x5102, 0x4223, 0x1326],
  [0x4223, 0x1326, 0x3304, 0x2315, 0x4113],
  [0x4113, 0x8337, 0x4223],
  [0x6201, 0x5102, 0x3304, 0x4223, 0x4113, 0x8337],
  [0x6201, 0x2315, 0x8337, 0x4223],
  [0x5102, 0x3304, 0x2315, 0x8337, 0x4223],
  [0x5102, 0x4113, 0x8337, 0x1326],
  [0x4113, 0x8337, 0x1326, 0x3304, 0x6201],
  [0x6201, 0x2315, 0x8337, 0x1326, 0x5102],
  [0x3304, 0x2315, 0x8337, 0x1326],
  [0x3304, 0x1146, 0x2245],
  [0x6201, 0x5102, 0x1146, 0x2245],
  [0x6201, 0x2315, 0x4113, 0x3304, 0x1146, 0x2245],
  [0x2315, 0x4113, 0x5102, 0x1146, 0x2245],
  [0x5102, 0x4223, 0x1326, 0x3304, 0x1146, 0x2245],
  [0x1146, 0x2245, 0x6201, 0x4223, 0x1326],
  [0x3304, 0x1146, 0x2245, 0x6201, 0x2315, 0x4113, 0x5102, 0x4223, 0x1326],
  [0x4223, 0x1326, 0x1146, 0x2245, 0x2315, 0x4113],
  [0x4223, 0x4113, 0x8337, 0x3304, 0x1146, 0x2245],
  [0x6201, 0x5102, 0x1146, 0x2245, 0x4223, 0x4113, 0x8337],
  [0x4223, 0x6201, 0x2315, 0x8337, 0x3304, 0x1146, 0x2245],
  [0x4223, 0x8337, 0x2315, 0x2245, 0x1146, 0x5102],
  [0x5102, 0x4113, 0x8337, 0x1326, 0x3304, 0x1146, 0x2245],
  [0x4113, 0x8337, 0x1326, 0x1146, 0x2245, 0x6201],
  [0x6201, 0x2315, 0x8337, 0x1326, 0x5102, 0x3304, 0x1146, 0x2245],
  [0x2245, 0x2315, 0x8337, 0x1326, 0x1146],
  [0x2315, 0x2245, 0x8157],
  [0x6201, 0x5102, 0x3304, 0x2315, 0x2245, 0x8157],
  [0x4113, 0x6201, 0x2245, 0x8157],
  [0x2245, 0x8157, 0x4113, 0x5102, 0x3304],
  [0x5102, 0x4223, 0x1326, 0x2315, 0x2245, 0x8157],
  [0x6201, 0x4223, 0x1326, 0x3304, 0x2315, 0x2245, 0x8157],
  [0x6201, 0x2245, 0x8157, 0x4113, 0x5102, 0x4223, 0x1326],
  [0x4223, 0x1326, 0x3304, 0x2245, 0x8157, 0x4113],
  [0x4223, 0x4113, 0x8337, 0x2315, 0x2245, 0x8157],
  [0x6201, 0x5102, 0x3304, 0x4223, 0x4113, 0x8337, 0x2315, 0x2245, 0x8157],
  [0x8337, 0x4223, 0x6201, 0x2245, 0x8157],
  [0x5102, 0x3304, 0x2245, 0x8157, 0x8337, 0x4223],
  [0x5102, 0x4113, 0x8337, 0x1326, 0x2315, 0x2245, 0x8157],
  [0x4113, 0x8337, 0x1326, 0x3304, 0x6201, 0x2315, 0x2245, 0x8157],
  [0x5102, 0x1326, 0x8337, 0x8157, 0x2245, 0x6201],
  [0x8157, 0x8337, 0x1326, 0x3304, 0x2245],
  [0x2315, 0x3304, 0x1146, 0x8157],
  [0x6201, 0x5102, 0x1146, 0x8157, 0x2315],
  [0x3304, 0x1146, 0x8157, 0x4113, 0x6201],
  [0x4113, 0x5102, 0x1146, 0x8157],
  [0x2315, 0x3304, 0x1146, 0x8157, 0x5102, 0x4223, 0x1326],
  [0x1326, 0x4223, 0x6201, 0x2315, 0x8157, 0x1146],
  [0x3304, 0x1146, 0x8157, 0x4113, 0x6201, 0x5102, 0x4223, 0x1326],
  [0x1326, 0x1146, 0x8157, 0x4113, 0x4223],
  [0x2315, 0x3304, 0x1146, 0x8157, 0x4223, 0x4113, 0x8337],
  [0x6201, 0x5102, 0x1146, 0x8157, 0x2315, 0x4223, 0x4113, 0x8337],
  [0x3304, 0x1146, 0x8157, 0x8337, 0x4223, 0x6201],
  [0x4223, 0x5102, 0x1146, 0x8157, 0x8337],
  [0x2315, 0x3304, 0x1146, 0x8157, 0x5102, 0x4113, 0x8337, 0x1326],
  [0x6201, 0x4113, 0x8337, 0x1326, 0x1146, 0x8157, 0x2315],
  [0x6201, 0x3304, 0x1146, 0x8157, 0x8337, 0x1326, 0x5102],
  [0x1326, 0x1146, 0x8157, 0x8337],
  [0x1326, 0x8267, 0x1146],
  [0x6201, 0x5102, 0x3304, 0x1326, 0x8267, 0x1146],
  [0x6201, 0x2315, 0x4113, 0x1326, 0x8267, 0x1146],
  [0x5102, 0x3304, 0x2315, 0x4113, 0x1326, 0x8267, 0x1146],
  [0x5102, 0x4223, 0x8267, 0x1146],
  [0x3304, 0x6201, 0x4223, 0x8267, 0x1146],
  [0x5102, 0x4223, 0x8267, 0x1146, 0x6201, 0x2315, 0x4113],
  [0x1146, 0x8267, 0x4223, 0x4113, 0x2315, 0x3304],
  [0x4113, 0x8337, 0x4223, 0x1326, 0x8267, 0x1146],
  [0x6201, 0x5102, 0x3304, 0x4223, 0x4113, 0x8337, 0x1326, 0x8267, 0x1146],
  [0x6201, 0x2315, 0x8337, 0x4223, 0x1326, 0x8267, 0x1146],
  [0x5102, 0x3304, 0x2315, 0x8337, 0x4223, 0x1326, 0x8267, 0x1146],
  [0x8267, 0x1146, 0x5102, 0x4113, 0x8337],
  [0x6201, 0x4113, 0x8337, 0x8267, 0x1146, 0x3304],
  [0x6201, 0x2315, 0x8337, 0x8267, 0x1146, 0x5102],
  [0x1146, 0x3304, 0x2315, 0x8337, 0x8267],
  [0x3304, 0x1326, 0x8267, 0x2245],
  [0x1326, 0x8267, 0x2245, 0x6201, 0x5102],
  [0x3304, 0x1326, 0x8267, 0x2245, 0x6201, 0x2315, 0x4113],
  [0x1326, 0x8267, 0x2245, 0x2315, 0x4113, 0x5102],
  [0x5102, 0x4223, 0x8267, 0x2245, 0x3304],
  [0x6201, 0x4223, 0x8267, 0x2245],
  [0x5102, 0x4223, 0x8267, 0x2245, 0x3304, 0x6201, 0x2315, 0x4113],
  [0x4113, 0x4223, 0x8267, 0x2245, 0x2315],
  [0x3304, 0x1326, 0x8267, 0x2245, 0x4223, 0x4113, 0x8337],
  [0x1326, 0x8267, 0x2245, 0x6201, 0x5102, 0x4223, 0x4113, 0x8337],
  [0x3304, 0x1326, 0x8267, 0x2245, 0x4223, 0x6201, 0x2315, 0x8337],
  [0x5102, 0x1326, 0x8267, 0x2245, 0x2315, 0x8337, 0x4223],
  [0x3304, 0x2245, 0x8267, 0x8337, 0x4113, 0x5102],
  [0x8337, 0x8267, 0x2245, 0x6201, 0x4113],
  [0x5102, 0x6201, 0x2315, 0x8337, 0x8267, 0x2245, 0x3304],
  [0x2315, 0x8337, 0x8267, 0x2245],
  [0x2315, 0x2245, 0x8157, 0x1326, 0x8267, 0x1146],
  [0x6201, 0x5102, 0x3304, 0x2315, 0x2245, 0x8157, 0x1326, 0x8267, 0x1146],
  [0x6201, 0x2245, 0x8157, 0x4113, 0x1326, 0x8267, 0x1146],
  [0x2245, 0x8157, 0x4113, 0x5102, 0x3304, 0x1326, 0x8267, 0x1146],
  [0x4223, 0x8267, 0x1146, 0x5102, 0x2315, 0x2245, 0x8157],
  [0x3304, 0x6201, 0x4223, 0x8267, 0x1146, 0x2315, 0x2245, 0x8157],
  [0x4223, 0x8267, 0x1146, 0x5102, 0x6201, 0x2245, 0x8157, 0x4113],
  [0x3304, 0x2245, 0x8157, 0x4113, 0x4223, 0x8267, 0x1146],
  [0x4223, 0x4113, 0x8337, 0x2315, 0x2245, 0x8157, 0x1326, 0x8267, 0x1146],
  [
    0x6201,
    0x5102,
    0x3304,
    0x4223,
    0x4113,
    0x8337,
    0x2315,
    0x2245,
    0x8157,
    0x1326,
    0x8267,
    0x1146
  ],
  [0x8337, 0x4223, 0x6201, 0x2245, 0x8157, 0x1326, 0x8267, 0x1146],
  [0x4223, 0x5102, 0x3304, 0x2245, 0x8157, 0x8337, 0x1326, 0x8267, 0x1146],
  [0x8267, 0x1146, 0x5102, 0x4113, 0x8337, 0x2315, 0x2245, 0x8157],
  [0x6201, 0x4113, 0x8337, 0x8267, 0x1146, 0x3304, 0x2315, 0x2245, 0x8157],
  [0x8337, 0x8267, 0x1146, 0x5102, 0x6201, 0x2245, 0x8157],
  [0x3304, 0x2245, 0x8157, 0x8337, 0x8267, 0x1146],
  [0x8157, 0x2315, 0x3304, 0x1326, 0x8267],
  [0x8267, 0x8157, 0x2315, 0x6201, 0x5102, 0x1326],
  [0x8267, 0x1326, 0x3304, 0x6201, 0x4113, 0x8157],
  [0x8267, 0x8157, 0x4113, 0x5102, 0x1326],
  [0x5102, 0x4223, 0x8267, 0x8157, 0x2315, 0x3304],
  [0x2315, 0x6201, 0x4223, 0x8267, 0x8157],
  [0x3304, 0x5102, 0x4223, 0x8267, 0x8157, 0x4113, 0x6201],
  [0x4113, 0x4223, 0x8267, 0x8157],
  [0x8157, 0x2315, 0x3304, 0x1326, 0x8267, 0x4223, 0x4113, 0x8337],
  [0x8157, 0x2315, 0x6201, 0x5102, 0x1326, 0x8267, 0x4223, 0x4113, 0x8337],
  [0x8157, 0x8337, 0x4223, 0x6201, 0x3304, 0x1326, 0x8267],
  [0x5102, 0x1326, 0x8267, 0x8157, 0x8337, 0x4223],
  [0x8267, 0x8157, 0x2315, 0x3304, 0x5102, 0x4113, 0x8337],
  [0x6201, 0x4113, 0x8337, 0x8267, 0x8157, 0x2315],
  [0x6201, 0x3304, 0x5102, 0x8337, 0x8267, 0x8157],
  [0x8337, 0x8267, 0x8157],
  [0x8337, 0x8157, 0x8267],
  [0x6201, 0x5102, 0x3304, 0x8337, 0x8157, 0x8267],
  [0x6201, 0x2315, 0x4113, 0x8337, 0x8157, 0x8267],
  [0x5102, 0x3304, 0x2315, 0x4113, 0x8337, 0x8157, 0x8267],
  [0x5102, 0x4223, 0x1326, 0x8337, 0x8157, 0x8267],
  [0x6201, 0x4223, 0x1326, 0x3304, 0x8337, 0x8157, 0x8267],
  [0x6201, 0x2315, 0x4113, 0x5102, 0x4223, 0x1326, 0x8337, 0x8157, 0x8267],
  [0x4223, 0x1326, 0x3304, 0x2315, 0x4113, 0x8337, 0x8157, 0x8267],
  [0x4113, 0x8157, 0x8267, 0x4223],
  [0x4223, 0x4113, 0x8157, 0x8267, 0x6201, 0x5102, 0x3304],
  [0x8157, 0x8267, 0x4223, 0x6201, 0x2315],
  [0x3304, 0x2315, 0x8157, 0x8267, 0x4223, 0x5102],
  [0x1326, 0x5102, 0x4113, 0x8157, 0x8267],
  [0x8157, 0x4113, 0x6201, 0x3304, 0x1326, 0x8267],
  [0x1326, 0x5102, 0x6201, 0x2315, 0x8157, 0x8267],
  [0x8267, 0x1326, 0x3304, 0x2315, 0x8157],
  [0x3304, 0x1146, 0x2245, 0x8337, 0x8157, 0x8267],
  [0x6201, 0x5102, 0x1146, 0x2245, 0x8337, 0x8157, 0x8267],
  [0x6201, 0x2315, 0x4113, 0x3304, 0x1146, 0x2245, 0x8337, 0x8157, 0x8267],
  [0x2315, 0x4113, 0x5102, 0x1146, 0x2245, 0x8337, 0x8157, 0x8267],
  [0x5102, 0x4223, 0x1326, 0x3304, 0x1146, 0x2245, 0x8337, 0x8157, 0x8267],
  [0x1146, 0x2245, 0x6201, 0x4223, 0x1326, 0x8337, 0x8157, 0x8267],
  [
    0x6201,
    0x2315,
    0x4113,
    0x5102,
    0x4223,
    0x1326,
    0x3304,
    0x1146,
    0x2245,
    0x8337,
    0x8157,
    0x8267
  ],
  [0x4113, 0x4223, 0x1326, 0x1146, 0x2245, 0x2315, 0x8337, 0x8157, 0x8267],
  [0x4223, 0x4113, 0x8157, 0x8267, 0x3304, 0x1146, 0x2245],
  [0x6201, 0x5102, 0x1146, 0x2245, 0x4223, 0x4113, 0x8157, 0x8267],
  [0x8157, 0x8267, 0x4223, 0x6201, 0x2315, 0x3304, 0x1146, 0x2245],
  [0x2315, 0x8157, 0x8267, 0x4223, 0x5102, 0x1146, 0x2245],
  [0x1326, 0x5102, 0x4113, 0x8157, 0x8267, 0x3304, 0x1146, 0x2245],
  [0x1326, 0x1146, 0x2245, 0x6201, 0x4113, 0x8157, 0x8267],
  [0x5102, 0x6201, 0x2315, 0x8157, 0x8267, 0x1326, 0x3304, 0x1146, 0x2245],
  [0x1326, 0x1146, 0x2245, 0x2315, 0x8157, 0x8267],
  [0x2315, 0x2245, 0x8267, 0x8337],
  [0x2315, 0x2245, 0x8267, 0x8337, 0x6201, 0x5102, 0x3304],
  [0x4113, 0x6201, 0x2245, 0x8267, 0x8337],
  [0x5102, 0x4113, 0x8337, 0x8267, 0x2245, 0x3304],
  [0x2315, 0x2245, 0x8267, 0x8337, 0x5102, 0x4223, 0x1326],
  [0x6201, 0x4223, 0x1326, 0x3304, 0x8337, 0x2315, 0x2245, 0x8267],
  [0x4113, 0x6201, 0x2245, 0x8267, 0x8337, 0x5102, 0x4223, 0x1326],
  [0x4113, 0x4223, 0x1326, 0x3304, 0x2245, 0x8267, 0x8337],
  [0x2315, 0x2245, 0x8267, 0x4223, 0x4113],
  [0x2315, 0x2245, 0x8267, 0x4223, 0x4113, 0x6201, 0x5102, 0x3304],
  [0x6201, 0x2245, 0x8267, 0x4223],
  [0x3304, 0x2245, 0x8267, 0x4223, 0x5102],
  [0x5102, 0x4113, 0x2315, 0x2245, 0x8267, 0x1326],
  [0x4113, 0x2315, 0x2245, 0x8267, 0x1326, 0x3304, 0x6201],
  [0x5102, 0x6201, 0x2245, 0x8267, 0x1326],
  [0x3304, 0x2245, 0x8267, 0x1326],
  [0x8267, 0x8337, 0x2315, 0x3304, 0x1146],
  [0x5102, 0x1146, 0x8267, 0x8337, 0x2315, 0x6201],
  [0x3304, 0x1146, 0x8267, 0x8337, 0x4113, 0x6201],
  [0x8337, 0x4113, 0x5102, 0x1146, 0x8267],
  [0x8267, 0x8337, 0x2315, 0x3304, 0x1146, 0x5102, 0x4223, 0x1326],
  [0x1146, 0x8267, 0x8337, 0x2315, 0x6201, 0x4223, 0x1326],
  [0x8267, 0x8337, 0x4113, 0x6201, 0x3304, 0x1146, 0x5102, 0x4223, 0x1326],
  [0x4113, 0x4223, 0x1326, 0x1146, 0x8267, 0x8337],
  [0x3304, 0x2315, 0x4113, 0x4223, 0x8267, 0x1146],
  [0x2315, 0x6201, 0x5102, 0x1146, 0x8267, 0x4223, 0x4113],
  [0x1146, 0x8267, 0x4223, 0x6201, 0x3304],
  [0x5102, 0x1146, 0x8267, 0x4223],
  [0x8267, 0x1326, 0x5102, 0x4113, 0x2315, 0x3304, 0x1146],
  [0x6201, 0x4113, 0x2315, 0x1326, 0x1146, 0x8267],
  [0x6201, 0x3304, 0x1146, 0x8267, 0x1326, 0x5102],
  [0x1326, 0x1146, 0x8267],
  [0x1326, 0x8337, 0x8157, 0x1146],
  [0x8337, 0x8157, 0x1146, 0x1326, 0x6201, 0x5102, 0x3304],
  [0x8337, 0x8157, 0x1146, 0x1326, 0x6201, 0x2315, 0x4113],
  [0x4113, 0x5102, 0x3304, 0x2315, 0x1326, 0x8337, 0x8157, 0x1146],
  [0x8337, 0x8157, 0x1146, 0x5102, 0x4223],
  [0x6201, 0x4223, 0x8337, 0x8157, 0x1146, 0x3304],
  [0x8337, 0x8157, 0x1146, 0x5102, 0x4223, 0x6201, 0x2315, 0x4113],
  [0x4223, 0x8337, 0x8157, 0x1146, 0x3304, 0x2315, 0x4113],
  [0x4223, 0x4113, 0x8157, 0x1146, 0x1326],
  [0x4223, 0x4113, 0x8157, 0x1146, 0x1326, 0x6201, 0x5102, 0x3304],
  [0x1146, 0x8157, 0x2315, 0x6201, 0x4223, 0x1326],
  [0x4223, 0x5102, 0x3304, 0x2315, 0x8157, 0x1146, 0x1326],
  [0x4113, 0x8157, 0x1146, 0x5102],
  [0x6201, 0x4113, 0x8157, 0x1146, 0x3304],
  [0x2315, 0x8157, 0x1146, 0x5102, 0x6201],
  [0x2315, 0x8157, 0x1146, 0x3304],
  [0x2245, 0x3304, 0x1326, 0x8337, 0x8157],
  [0x6201, 0x2245, 0x8157, 0x8337, 0x1326, 0x5102],
  [0x2245, 0x3304, 0x1326, 0x8337, 0x8157, 0x6201, 0x2315, 0x4113],
  [0x2245, 0x2315, 0x4113, 0x5102, 0x1326, 0x8337, 0x8157],
  [0x4223, 0x8337, 0x8157, 0x2245, 0x3304, 0x5102],
  [0x8157, 0x2245, 0x6201, 0x4223, 0x8337],
  [0x2245, 0x3304, 0x5102, 0x4223, 0x8337, 0x8157, 0x4113, 0x6201, 0x2315],
  [0x4223, 0x8337, 0x8157, 0x2245, 0x2315, 0x4113],
  [0x4113, 0x8157, 0x2245, 0x3304, 0x1326, 0x4223],
  [0x1326, 0x4223, 0x4113, 0x8157, 0x2245, 0x6201, 0x5102],
  [0x8157, 0x2245, 0x3304, 0x1326, 0x4223, 0x6201, 0x2315],
  [0x5102, 0x1326, 0x4223, 0x2315, 0x8157, 0x2245],
  [0x3304, 0x5102, 0x4113, 0x8157, 0x2245],
  [0x4113, 0x8157, 0x2245, 0x6201],
  [0x5102, 0x6201, 0x2315, 0x8157, 0x2245, 0x3304],
  [0x2315, 0x8157, 0x2245],
  [0x1146, 0x1326, 0x8337, 0x2315, 0x2245],
  [0x1146, 0x1326, 0x8337, 0x2315, 0x2245, 0x6201, 0x5102, 0x3304],
  [0x6201, 0x2245, 0x1146, 0x1326, 0x8337, 0x4113],
  [0x2245, 0x1146, 0x1326, 0x8337, 0x4113, 0x5102, 0x3304],
  [0x5102, 0x1146, 0x2245, 0x2315, 0x8337, 0x4223],
  [0x1146, 0x3304, 0x6201, 0x4223, 0x8337, 0x2315, 0x2245],
  [0x8337, 0x4113, 0x6201, 0x2245, 0x1146, 0x5102, 0x4223],
  [0x4223, 0x8337, 0x4113, 0x3304, 0x2245, 0x1146],
  [0x4113, 0x2315, 0x2245, 0x1146, 0x1326, 0x4223],
  [0x1146, 0x1326, 0x4223, 0x4113, 0x2315, 0x2245, 0x6201, 0x5102, 0x3304],
  [0x1326, 0x4223, 0x6201, 0x2245, 0x1146],
  [0x4223, 0x5102, 0x3304, 0x2245, 0x1146, 0x1326],
  [0x2245, 0x1146, 0x5102, 0x4113, 0x2315],
  [0x4113, 0x2315, 0x2245, 0x1146, 0x3304, 0x6201],
  [0x6201, 0x2245, 0x1146, 0x5102],
  [0x3304, 0x2245, 0x1146],
  [0x3304, 0x1326, 0x8337, 0x2315],
  [0x5102, 0x1326, 0x8337, 0x2315, 0x6201],
  [0x6201, 0x3304, 0x1326, 0x8337, 0x4113],
  [0x5102, 0x1326, 0x8337, 0x4113],
  [0x4223, 0x8337, 0x2315, 0x3304, 0x5102],
  [0x6201, 0x4223, 0x8337, 0x2315],
  [0x3304, 0x5102, 0x4223, 0x8337, 0x4113, 0x6201],
  [0x4113, 0x4223, 0x8337],
  [0x4113, 0x2315, 0x3304, 0x1326, 0x4223],
  [0x1326, 0x4223, 0x4113, 0x2315, 0x6201, 0x5102],
  [0x3304, 0x1326, 0x4223, 0x6201],
  [0x5102, 0x1326, 0x4223],
  [0x5102, 0x4113, 0x2315, 0x3304],
  [0x6201, 0x4113, 0x2315],
  [0x6201, 0x3304, 0x5102],
  [],
];

const int vertexDataRowLength = 12;

final Uint16List vertexData =
    new Uint16List.fromList(_Flatten(_vertexDataOriginal, vertexDataRowLength));

void SanityCheckTables() {
  assert(256 == hotCornerMap.length);
  assert(256 == vertexData.length);
  for (int i = 0; i < hotCornerMap.length; i++) {
    final int cls = hotCornerMap[i];
    assert(cls < _configClassOriginal.length);
    List<int> cfg = _configClassOriginal[cls];
    assert(cfg[1] * 3 + 2 == cfg.length);

    assert(cfg[0] == _vertexDataOriginal[i].length);
    for (int j = 2; j < cfg.length; j++) {
      assert(cfg[j] < _vertexDataOriginal[i].length);
    }
  }
}
