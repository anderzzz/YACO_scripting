'''Class to read image data from Excel sheets and associate them with a line item. The methods
build on a specific layout of the Excel sheets wherein they have ordered images

'''
import zipfile
import pandas as pd
import xml.etree.ElementTree as ET
from PIL import Image

from enum import Enum

class WBArchive(Enum):
    '''An Enum with data to navigate Excel XML rawfiles. This abstracts details of how the Excel sheet is
    constructed from the methods.

    '''
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
    '''Parse and associate image content in Excel sheet

    The initialization of the class generates folder `xl` containing unzipped data from the Excel sheet. Under
    normal execution, these files and folders should be untouched by the user. After execution the folders can
    be manually deleted.

    Args:
        rawdata_path (str): Path to Excel sheet to process
        image_suffix (list, optional): Suffixes for the images to consider in the Excel sheet

    Attributes:
        sheetname2img (dict): In the rawdata, the association between data sheet and the image. Note that the same image
            can appear in multiple places.
        n_imgs (int): The total number of images found in the entire Excel sheet.

    '''
    def __init__(self, rawdata_path, image_suffix=['png','jpg']):

        self.rawdata_path = rawdata_path
        self.image_suffix = image_suffix
        self.association = None
        self.excel_archive = WBArchive
        self.excel_archive_namespaces = {'ns1' : self.excel_archive.XMLNameSpace1.value,
                                         'ns2' : self.excel_archive.XMLNameSpace2.value}

        # Unzip all image files in the Excel workbook
        self.img_paths = []
        for embedded_file in zipfile.ZipFile(self.rawdata_path).namelist():
            for img_type in self.image_suffix:
                if '.{}'.format(img_type) in embedded_file:
                    self.img_paths.append(zipfile.ZipFile(self.rawdata_path).extract(embedded_file))

        self.sheetname2img, self.sheetname2sheetid, self.sheetid2drawing, self.drawing2img = self._associate_sheets_imgs()

    def _associate_sheets_imgs(self):
        '''By parsing the XML raw data, an association between sheet name and rawdata image name is inferred.
        The association is found via intermediate associations. This function builds on details of how Excel stores
        its data in XML files. If different Excel versions are used, this function can break.

        '''
        # First associate a sheet name with a sheet ID. The association is found in the Workbook rawdata
        sheetname2id = {}
        with zipfile.ZipFile(self.rawdata_path) as zipper:
            tree = self._get_xml(zipper, self.excel_archive.WorkBookPath.value)
            for xx in tree.findall('ns1:sheets', self.excel_archive_namespaces):
                for sheet in xx:
                    sheetname2id[sheet.attrib[self.excel_archive.WorkBookSheetName.value]] = \
                        int(sheet.attrib[self.excel_archive.WorkBookSheetId.value])

        # Second associate a sheet ID with a drawing ID. A drawing is a collection of graphics that exists in a sheet
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

        # Third associate drawing ID with one image ID or a plurality of image IDs
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
        '''Return the rawdata images for a given sheet

        '''
        try:
            imgs = self.sheetname2img[sheet_name]
        except KeyError:
            raise KeyError('Did not find sheet with name {}'.format(sheet_name))

        return imgs

    def associate_imgs(self, sheet_name):
        '''Create the association between raw data images and experimental identifier.

        The method is not guaranteed to work and builds on reasonable assumptions of how the Excel sheet was constructed.
        If association is unsuccessful, error messages are printed, but no exception is raised.

        Args:
            sheet_name (str): The name of the sheet for which to create an association

        Returns:
            association (dict): Association between experimental identifier and the rawdata image name

        '''
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
        '''Save rawdata images that have been associated with respective experimental identifiers. This can only be
        run after the `associate_imgs` has been run

        Args:
            out_name (str): Path to an existing directory to which the renamed image files are written

        '''
        if self.association is None:
            raise RuntimeError('No image association found. Did you run `associate_imgs` prior?')

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


content = SheetContent('./data/29762_204A_Octet_2020Sep03.xlsx')
print (content.sheetname2img)
content.associate_imgs('DataTable')
content.save_assocation('dummy')
