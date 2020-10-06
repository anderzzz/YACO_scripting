'''Excel sheet reader

'''
import zipfile
import pandas as pd
import xml.etree.ElementTree as ET
from PIL import Image

from enum import Enum

class WBArchive(Enum):
    WorkBookPath = 'xl/workbook.xml'
    WorkBookSheetName = 'name'
    WorkBookSheetId = 'sheetId'
    WorkSheetRelsPathRoot = 'xl/worksheets/_rels/'
    WorkSheetName = 'sheet'
    WorkSheetRelsSuffix = '.xml.rels'
    WorkSheetRelsDrawing = 'Target'
    DrawingRelsPathRoot = 'xl/drawings/_rels/'
    DrawingRelsSuffix = '.rels'
    XMLNameSpace1 = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    XMLNameSpace2 = 'http://schemas.openxmlformats.org/package/2006/relationships'

class SheetContent(object):

    def __init__(self, rawdata_path, sheet_name_selection=None, image_suffix=['png','jpg']):

        self.rawdata_path = rawdata_path
        self.image_suffix = image_suffix
        self.excel_archive = WBArchive
        self.excel_archive_namespaces = {'ns1' : self.excel_archive.XMLNameSpace1.value,
                                         'ns2' : self.excel_archive.XMLNameSpace2.value}

        self.img_paths = []
        for embedded_file in zipfile.ZipFile(self.rawdata_path).namelist():
            for img_type in self.image_suffix:
                if '.{}'.format(img_type) in embedded_file:
                    self.img_paths.append(zipfile.ZipFile(self.rawdata_path).extract(embedded_file))

        self.sheetname2img, self.sheetname2sheetid, self.sheetid2drawing, self.drawing2img = self._associate_sheets_imgs()

    def _associate_sheets_imgs(self):

        # First associate a sheet name with a sheet ID
        sheetname2id = {}
        with zipfile.ZipFile(self.rawdata_path) as zipper:
            tree = self._get_xml(zipper, self.excel_archive.WorkBookPath.value)
            for xx in tree.findall('ns1:sheets', self.excel_archive_namespaces):
                for sheet in xx:
                    sheetname2id[sheet.attrib[self.excel_archive.WorkBookSheetName.value]] = \
                        int(sheet.attrib[self.excel_archive.WorkBookSheetId.value])

        # Second associate a sheet ID with a drawing ID
        sheetid2drawingid = {}
        with zipfile.ZipFile(self.rawdata_path) as zipper:
            for nsheet, _ in enumerate(sheetname2id):
                tree = self._get_xml(zipper,
                                     '{}{}{}{}'.format(self.excel_archive.WorkSheetRelsPathRoot.value,
                                                       self.excel_archive.WorkSheetName.value,
                                                       nsheet + 1,
                                                       self.excel_archive.WorkSheetRelsSuffix.value))
                for xx in tree.findall('ns2:Relationship', self.excel_archive_namespaces):
                    drawing_str = xx.attrib[self.excel_archive.WorkSheetRelsDrawing.value]
                    sheetid2drawingid[nsheet + 1] = drawing_str.split('/')[-1]

        # Third associate drawing ID with a plurality of image IDs
        drawingid2imageids = {}
        with zipfile.ZipFile(self.rawdata_path) as zipper:
            for drawing in sheetid2drawingid.values():
                tree = self._get_xml(zipper,
                                     '{}{}{}'.format(self.excel_archive.DrawingRelsPathRoot.value,
                                                     drawing,
                                                     self.excel_archive.DrawingRelsSuffix.value))
                imgs= []
                for xx in tree.findall('ns2:Relationship', self.excel_archive_namespaces):
                    img_str = xx.attrib[self.excel_archive.WorkSheetRelsDrawing.value]
                    imgs.append(img_str.split('/')[-1])
                drawingid2imageids[drawing] = imgs

        # Construct the sheet name to image ID map
        ret = {}
        for sheetname, id in sheetname2id.items():
            ret[sheetname] = drawingid2imageids[sheetid2drawingid[id]]

        return ret, sheetname2id, sheetid2drawingid, drawingid2imageids

    def _get_xml(self, zipper, path):
        with zipper.open(path) as fin:
            d1 = fin.read()
            return ET.fromstring(d1.decode('utf-8'))

    @property
    def n_imgs(self):
        return len(self.img_paths)

    def get_image_content(self, sheet_name):

        try:
            imgs = self.sheetname2img[sheet_name]
        except KeyError:
            raise KeyError('Did not find sheet with name {}'.format(sheet_name))

        return imgs

    def associate_imgs(self, sheet_name):

        def ordered_unique(seq):
            seen = set()
            seen_add = seen.add
            return [x for x in seq if not(x in seen or seen_add(x))]

        df = pd.read_excel(self.rawdata_path, sheet_name)
        vars = ordered_unique(df['Group Key'].tolist())
        imgs = self.get_image_content(sheet_name)

        self.association = {}
        if len(vars) == len(imgs):
            for var, img in zip(vars, imgs):
                self.association[var] = img

        else:
            print ('Association failed. Images ({} in total) != Group Keys ({} in total).'.format(len(imgs), len(vars)))
            print ('Group keys: {}'.format(vars))

        return self.association

    def save_assocation(self, out_name):

        print (self.association)
        print (self.img_paths)
        for run, img_name in self.association.items():
            new_img_name = '{}.{}'.format(run, img_name.split('.')[-1])
            save_path = '{}/{}'.format(out_name, new_img_name)

            img_source = [x for x in self.img_paths if img_name in x]
            if len(img_source) > 1:
                raise RuntimeError('Found multiple images that match the name {}'.format(img_name))
            elif len(img_source) == 0:
                raise RuntimeError('Found no image that matches the name {}'.format(img_name))

            img_data = Image.open(img_source[0])
            img_data = img_data.save(save_path)


ss = SheetContent('./data/29762_204A_Octet_2020Sep03.xlsx', sheet_name_selection='DataTable')
print (ss.n_imgs)
print (ss.get_image_content('DataTable'))
ss.associate_imgs('DataTable')
ss.save_assocation('dummy')
