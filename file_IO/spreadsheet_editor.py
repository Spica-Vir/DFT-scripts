#!/usr/bin/env python
# coding: utf-8

# A Spreadsheet_editor class to read/write common spreadsheet files
# Requirement: pandas package, imported as pd
# By Spica. Vir., ICL, Apr. 22, 22. spica.h.zhou@gmail.com

import pandas as pd


class Spreadsheet_editor:
    """
    A Spreadsheel_editor class:
        1. Open and read the common spreadsheet files
        2. Read user-specified rows and columns, export data as a list variable
    By Spica. Vir., ICL, Apr. 22, 22. spica.h.zhou@gmail.com
    """

    def __init__(self, file, sheet=0):
        """
        initialization & check
        choose_function: dict, save the available file format - function pairs
        file: str, file name with extension and, if necessary, full directory
        extension: str, file extension
        sheet: int/str, name or index of a sheet in a spreadsheet file
        """
        self.choose_function = {
            'xls': pd.read_excel,
            'xlsx': pd.read_excel,
            'csv': pd.read_csv,
        }
        extension = file.split('.')
        self.file = file
        self.extension = extension[-1]
        self.sheet = sheet

        if not self.choose_function.__contains__(self.extension):
            print('Error: File extension not supported.')
            return

    def read_data(self, row=[], col=[], data_format='list'):
        """
        Read the data of the sheet or of the specified rows/columns
        row: list, indices of rows to read
        col: list, indices of columns to read
        data_format: str, 'list' or any other string, the format of output data
                    available choices: list / pandas dataframe

        IMPORTANT NOTES:
        1. Indices of rows and columns begin from 1, i.e., consistent with
           indices in the original file.
        2. Any header/index in the sheet will be recognized as values.
           Excluding them when setting 'row' and 'col', or after exporting data
           is needed.
        """
        if self.sheet:
            complete_data = self.choose_function[self.extension](
                self.file, sheet_name=self.sheet, header=None)
        else:
            complete_data = self.choose_function[self.extension](
                self.file, header=None)

        data_out = complete_data
        if row:
            row = [r - 1 for r in row]
            data_out = data_out.iloc[row, :]

        if col:
            col = [c - 1 for c in col]
            data_out = data_out.iloc[:, col]

        if data_format == 'list':
            data_out = data_out.values.tolist()

        return data_out
