import re
import tempfile
import zipfile
from io import BytesIO
from pathlib import Path

import FreeCAD


def get_material_and_movement(label):
    start = label.find("(")
    end = label.find(")")
    if (start == -1) or (end == -1):
        return None
    # print(label,start,end)
    material_and_movement = label[start + 1:end]
    comma_pos = material_and_movement.find(',')
    if comma_pos > 0:
        has_movement = True
        name = material_and_movement[:comma_pos]
    else:
        has_movement = False
        name = material_and_movement
    return name, has_movement


def test_freecad_file_ok(filename):
    try:
        FreeCAD.openDocument(filename)
        doc = FreeCAD.ActiveDocument
        FreeCAD.closeDocument(doc.Name)
        return True
    except:
        return False


def test_freecad_file_bytes_ok(file_bytes):
    with tempfile.NamedTemporaryFile('wb') as tf:
        tf.write(file_bytes)
        filename = tf.name
        return test_freecad_file_ok(filename)


def get_freecad_filenames_from_zipobject(zipobject):
    freecad_filenames = []
    for fn in zipobject.filelist:
        fb = zipobject.read(fn)
        if test_freecad_file_bytes_ok(fb):
            freecad_filenames.append(fn)
    return freecad_filenames


def get_freecad_filenames_from_zip(filename):
    freecad_filenames = []
    try:
        with zipfile.ZipFile(filename) as zf:
            for fn in zf.filelist:
                fb = zf.read(fn)
                if test_freecad_file_bytes_ok(fb):
                    freecad_filenames.append(fn)
    except:
        pass
    return freecad_filenames


def test_zip_of_freecad_files_ok(filename):
    return len(get_freecad_filenames_from_zip(filename)) > 0


def test_zip_of_freecad_files_bytes_ok(file_bytes):
    with tempfile.NamedTemporaryFile('wb') as tf:
        tf.write(file_bytes)
        filename = tf.name
        return test_zip_of_freecad_files_ok(filename)


def materials_for_zipfile_bytes(zipfile_bytes):
    materials_labels = set()
    zfobject = zipfile.ZipFile(BytesIO(zipfile_bytes))
    freecad_filenames = get_freecad_filenames_from_zipobject(zfobject)
    for fn in freecad_filenames:
        file_bytes = zfobject.read(fn)
        new_mats = materials_for_file_bytes(file_bytes)
        materials_labels |= new_mats
    return materials_labels


def materials_for_file_bytes(file_bytes):
    with tempfile.NamedTemporaryFile('wb') as tf:
        tf.write(file_bytes)
        filename = tf.name
        mats = materials_for_filename(filename)
    return mats


def materials_for_filename(filename):
    material_labels = set()
    FreeCAD.openDocument(str(filename))
    doc = FreeCAD.ActiveDocument
    objects = doc.Objects

    for obj in objects:
        # noinspection PyNoneFunctionAssignment
        label = obj.Label
        mat_mov = get_material_and_movement(label)
        if mat_mov is None:
            continue
        name, _ = mat_mov
        material_labels.add(name)

    extra_params = parameters_for_freecad_doc(doc)
    if 'vacuum_material' in extra_params:
        material_labels.add(extra_params['vacuum_material'][0])

    FreeCAD.closeDocument(doc.Name)
    return material_labels


re_digits = re.compile('\d+')
re_letters = re.compile('[A-Z]+')


def get_non_empty_cells(sheet):
    import xml.etree.ElementTree as ET
    root = ET.fromstring(sheet.Content)
    return [c.attrib['address'] for c in root.iter('Cell')]


def get_col_row(cell):
    """given AAANNN, return NNN"""
    return re_letters.search(cell).group(), re_digits.search(cell).group()


def get_document(doc: FreeCAD.Document | str | Path):
    if type(doc) is Path:
        doc = str(doc)
    if type(doc) is str:
        doc = FreeCAD.openDocument(doc)
    return doc


def parameters_for_freecad_doc(doc: FreeCAD.Document) -> dict:
    try:
        sheet = doc.Spreadsheet
    except AttributeError:
        return {}
    cell_addresses = get_non_empty_cells(sheet)
    col_row_addresses = list(map(get_col_row, cell_addresses))
    # rows = set([cr[0] for cr in col_row_addresses])
    parameters = {}
    # get identifiers
    for (c, r) in col_row_addresses:
        if c != 'A':
            continue
        identifier = sheet.get(c+r)
        parameters[identifier] = []
    # get values
    for (c, r) in col_row_addresses:
        if c == 'A':
            continue
        identifier = sheet.get('A'+r)
        value = sheet.get(c+r)
        parameters[identifier].append(value)
    return parameters


def get_aperture_pv(doc):
    try:
        return parameters_for_freecad_doc(doc)['aperture_pv'][0]
    except:
        return None


def get_aperture_th(doc):
    try:
        return parameters_for_freecad_doc(doc)['aperture_th'][0]
    except:
        return None

