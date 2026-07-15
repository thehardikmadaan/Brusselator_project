import xml.etree.ElementTree as ET
import base64
import zlib
import struct

def read_vts_xml(filepath):
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    appended_data = root.find('AppendedData').text
    if appended_data.startswith('_'):
        appended_data = appended_data[1:]
    
    appended_data = ''.join(appended_data.split())
    
    arrays = []
    for da in root.findall('.//DataArray'):
        if da.get('format') == 'appended':
            arrays.append((da.get('Name'), int(da.get('offset'))))
            
    arrays.sort(key=lambda x: x[1])
    
    def decode_array(start_idx, end_idx=None):
        b64_str = appended_data[start_idx:end_idx] if end_idx else appended_data[start_idx:]
        raw_data = base64.b64decode(b64_str)
        
        num_blocks = struct.unpack_from('<I', raw_data, 0)[0]
        header_fmt = '<' + str(3 + num_blocks) + 'I'
        header_size = struct.calcsize(header_fmt)
        headers = struct.unpack_from(header_fmt, raw_data, 0)
        
        comp_sizes = headers[3:]
        uncomp_data = bytearray()
        current_offset = header_size
        
        for comp_size in comp_sizes:
            comp_block = raw_data[current_offset : current_offset + comp_size]
            uncomp_data.extend(zlib.decompress(comp_block))
            current_offset += comp_size
            
        return len(uncomp_data)

    c1_idx = next(i for i, (name, _) in enumerate(arrays) if name == 'C1')
    c2_idx = next(i for i, (name, _) in enumerate(arrays) if name == 'C2')
    
    c1_end = arrays[c1_idx+1][1] if c1_idx + 1 < len(arrays) else None
    c2_end = arrays[c2_idx+1][1] if c2_idx + 1 < len(arrays) else None
    
    c1 = decode_array(arrays[c1_idx][1], c1_end)
    c2 = decode_array(arrays[c2_idx][1], c2_end)
    
    return c1, c2

print(read_vts_xml("tests/final_output_30s (1).vts"))
