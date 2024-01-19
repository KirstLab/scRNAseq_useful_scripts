#!/usr/bin/env python
# coding: utf-8

import loompy

files = ["path/file1.loom",
         "path/file2.loom",
         "path/file3.loom",
         "path/file4.loom",
         "path/file5.loom"]

loompy.combine(files, 'path/your_output)file.loom')