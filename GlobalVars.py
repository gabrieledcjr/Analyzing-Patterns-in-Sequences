#!/usr/bin/python -tt

# Global flag variable for debugging mode
DEBUG = False

# Character shifting to left to compare for similarity
SHIFT = 5

HTML_HEAD_TAG = """
<head>
  <style type="text/css">
    .A {
      color: hotpink;
      font-family: monospace;
      font-size: 10px;
    }
    .B {
      color: olivedrab;
      font-family: monospace;
      font-size: 10px;
    }
    .C {
      color: orchid;
      font-family: monospace;
      font-size: 10px;
    }
    .D {
      color: mediumvioletred;
      font-family: monospace;
      font-size: 10px;
    }
    .E {
      color: lightseagreen;
      font-family: monospace;
      font-size: 10px;
    }
    .F {
      color: indigo;
      font-family: monospace;
      font-size: 10px;
    }
    .G {
      color: lightcoral;
      font-family: monospace;
      font-size: 10px;
    }
    .H {
      color: red;
      font-family: monospace;
      font-size: 10px;
    }
    .I {
      color: orangered;
      font-family: monospace;
      font-size: 10px;
    }
    .K {
      color: steelblue;
      font-family: monospace;
      font-size: 10px;
    }
    .L {
      color: red;
      font-family: monospace;
      font-size: 10px;
    }
    .M {
      color: green;
      font-family: monospace;
      font-size: 10px;
    }
    .N {
      color: blue;
      font-family: monospace;
      font-size: 10px;
    }
    .P {
      color: darkgoldenrod;
      font-family: monospace;
      font-size: 10px;
    }
    .Q {
      color: maroon;
      font-family: monospace;
      font-size: 10px;
    }
    .R {
      color: midnightblue;
      font-family: monospace;
      font-size: 10px;
    }
    .S {
      color: darkkhaki;
      font-family: monospace;
      font-size: 10px;
    }
    .T {
      color: olive;
      font-family: monospace;
      font-size: 10px;
    }
    .V {
      color: magenta;
      font-family: monospace;
      font-size: 10px;
    }
    .W {
      color: rosybrown;
      font-family: monospace;
      font-size: 10px;
    }
    .X {
      color: tomato;
      font-family: monospace;
      font-size: 10px;
    }
    .Y {
      color: gold;
      font-family: monospace;
      font-size: 10px;
    }
    .Z {
      color: royalblue;
      font-family: monospace;
      font-size: 10px;
    }
    .* {
      color: slategray;
      font-family: monospace;
      font-size: 10px;
    }
    .TEXT {
      color: black;
      font-family: monospace;
      font-size: 10px;
    }
    .SPACE {
      color: white;
      font-family: monospace;
      font-size: 10px;
    }
  </style>
<head>
"""
