#!/usr/bin/env python
#-*- coding:utf-8 -*-

import xlrd
import xlwt
from datetime import date,datetime

def read_excle():
    workbook = xlrd.open_workbook('Methylated_RNA_sites.mRNA.xlsx')
    # print(workbook.sheet_names())
    sheet2_name = workbook.sheet_names()[1]
    print(sheet2_name)

    sheet1 = workbook.sheets()[0]
    sheet2 = workbook.sheet_by_index(1)
    sheet3 = workbook.sheet_by_name('AIR3')
    print(sheet2,sheet3)

    print(sheet2.name,sheet2.nrows,sheet2.ncols)

    rows = sheet2.row_values(24)
    cols = sheet2.col_values(2)

    # print(rows,cols)
    print(rows)

    print(sheet2.cell(24,3),sheet2.cell(24,3).value)
    print(sheet2.cell_value(23,6))
    print(sheet2.row(23)[2],sheet2.row(23)[2].value)

    print(sheet2.cell(23,6).ctype)

read_excle()



def set_style(name, height, bold=False):
    style = xlwt.XFStyle()  # 初始化样式

    font = xlwt.Font()  # 为样式创建字体
    font.name = name  # 'Times New Roman'
    font.bold = bold
    # f.underline= Font.UNDERLINE_DOUBLE
    font.color_index = 4
    font.height = height

    # borders= xlwt.Borders()
    # borders.left= 6
    # borders.right= 6
    # borders.top= 6
    # borders.bottom= 6

    style.font = font
    # style.borders = borders

    return style




# 写excel
def write_excel():
    f = xlwt.Workbook()  # 创建工作簿

    '''
    创建第一个sheet:
        sheet1
    '''
    sheet1 = f.add_sheet(u'sheet1', cell_overwrite_ok=True)  # 创建sheet
    row0 = [u'业务', u'状态', u'北京', u'上海', u'广州', u'深圳', u'状态小计', u'合计']
    column0 = [u'机票', u'船票', u'火车票', u'汽车票', u'其它']
    status = [u'预订', u'出票', u'退票', u'业务小计']

    # 生成第一行
    for i in range(0, len(row0)):
        sheet1.write(0, i, row0[i], set_style('Times New Roman', 220, True))

    # 生成第一列和最后一列(合并4行)
    i, j = 1, 0
    while i < 4 * len(column0) and j < len(column0):
        sheet1.write_merge(i, i + 3, 0, 0, column0[j], set_style('Arial', 220, True))  # 第一列
        sheet1.write_merge(i, i + 3, 7, 7)  # 最后一列"合计"
        i += 4
        j += 1

    sheet1.write_merge(21, 21, 0, 1, u'合计', set_style('Times New Roman', 220, True))

    # 生成第二列
    i = 0
    while i < 4 * len(column0):
        for j in range(0, len(status)):
            sheet1.write(j + i + 1, 1, status[j])
        i += 4

    '''
    创建第二个sheet:
        sheet2
    '''
    sheet2 = f.add_sheet(u'sheet2', cell_overwrite_ok=True)  # 创建sheet
    row0 = [u'姓名', u'年龄', u'出生日期', u'爱好', u'关系']
    column0 = [u'小杰', u'小胖', u'小明', u'大神', u'大仙', u'小敏', u'无名']

    # 生成第一行
    for i in range(0, len(row0)):
        sheet2.write(0, i, row0[i], set_style('Times New Roman', 220, True))

    # 生成第一列
    for i in range(0, len(column0)):
        sheet2.write(i + 1, 0, column0[i], set_style('Times New Roman', 220))

    sheet2.write(1, 2, '1991/11/11')
    sheet2.write_merge(7, 7, 2, 4, u'暂无')  # 合并列单元格
    sheet2.write_merge(1, 2, 4, 4, u'好朋友')  # 合并行单元格

    # 插入图片
    # sheet2.insert_bitmap('1.bmp', 10, 2)

    # 添加超链接
    n = "HYPERLINK"
    sheet2.write_merge(9, 9, 2, 8, xlwt.Formula(n + '("http://www.cnblogs.com/zhoujie";"jzhou\'s blog")'),
                       set_style('Arial', 300, True))
    sheet2.write_merge(10, 10, 2, 8, xlwt.Formula(n + '("mailto:zhoujie0111@126.com";"contact me")'),
                       set_style('Arial', 300, True))

    f.save('demo1.xlsx')  # 保存文件


if __name__ == '__main__':
   # read_excel()
    write_excel()
