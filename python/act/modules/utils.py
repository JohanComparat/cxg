import pdb

def suffix(in_str):
    """ return suffix of (file) in_str"""

    from os.path import splitext
    return splitext(in_str)[1]

def dm_free_metadata_reader(xml_file, tags_chain):
    """
    General XML reader: 
    return the value from tags.chain in xml_file
    """
    
    import xml.etree.ElementTree as ET
    try:             
        with open(xml_file, 'rt') as fi:
            tree = ET.parse(fi)
    except UnicodeDecodeError:
        pdb.set_trace()
        print('xml_file is impossible to parse')
        tree = None
    current = tree.getroot()
    for tag in tags_chain.split("."):
        current = current.find(tag)
        
    return current.text
