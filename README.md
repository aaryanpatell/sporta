SPORTA is developed as a spot detection and screening software for X-ray diffraction data for HDF5. The program was originally developed for MX Industry at the Canadian Light Source (CLS), but has been open-sourced for extended public use. The program also includes a python wrapper for data into excel for analysis or reporting.

## INSTALL

Build dependencies in 'src' directory through the 'Makefile':
```bash
cd src; make
```
Build 'main.c' in main directory using 'Makefile':
```bash
cd ..; make
```

## COMMAND LINE INTERFACE (CLI)

The program accepts the following command-line options:

- `-h`: Display the help message.
- `-f <file>`: Specify the input file.
- `-o`: Enable output.
- `-i <index>`: Set the index to a specified value.
- `-d <data>`: Set the data to a specified value.
- `-r <dir>`: Enable traversing on path `<dir>`.
- `-t <value>`: Set the threshold scale factor to a specified value.
- `-s <dir>`: Set the search directory to a specified value.
- `-e <value>`: Set the ring exclusion proximity radius.
- `-p`: Enable histogram through GNUPLOT.
- `-g <sigma>`: Set gaussian filtering sigma value.
- `-n`: Display the notes, tips, and comments.

## USAGE

Example usage and arugment passing:

```bash
./sporta -f input_file.h5 -o
```

## LICENSE

Copyright (c) 2023, Aaryan Patel

MIT LICENSE

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Acknowledgments: Canadian Light Source, Industrial Science Group (CLS), Denis Spasyuk (CLS)
