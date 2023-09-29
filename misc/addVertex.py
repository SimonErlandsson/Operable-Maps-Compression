#WKT
 wkt = bin.decode('utf-8') # Decompressing data
        point_idx = 0
        prev = ''
        for c_idx, char in enumerate(wkt):
            if char == ',' and prev != ')':
                if insert_idx == point_idx:
                    insert_string = f', {pos[0]} {pos[1]}'
                    wkt = wkt[:c_idx] + insert_string + wkt[c_idx:]
                    break
                point_idx += 1
            prev = char
        bin = bytes(wkt, 'utf-8')


#WKTCOMP
        wkt = gzip.decompress(bin).decode('utf-8') # Decompressing data
        point_idx = 0
        prev = ''
        for c_idx, char in enumerate(wkt):
            if char == ',' and prev != ')':
                if insert_idx == point_idx:
                    insert_string = f', {pos[0]} {pos[1]}'
                    wkt = wkt[:c_idx] + insert_string + wkt[c_idx:]
                    break
                point_idx += 1
            prev = char
        bin = gzip.compress(bytes(wkt, 'utf-8')) 

# WKB, MetaWKB, WkbComp
        _, geometry = self.decompress(bin)
        wkt = shapely.to_wkt(geometry, rounding_precision=-1)
        point_idx = 0
        prev = ''
        for c_idx, char in enumerate(wkt):
            if char == ',' and prev != ')':
                if insert_idx == point_idx:
                    insert_string = f', {pos[0]} {pos[1]}'
                    wkt = wkt[:c_idx] + insert_string + wkt[c_idx:]
                    break
                point_idx += 1
            prev = char
        geometry = shapely.wkt.loads(wkt)
        _, bin = self.compress(geometry)

# MetaWKT
data = bin.split(b'\n', 1)[1] # Split at newline byte
        wkt = gzip.decompress(data).decode('utf-8')
        point_idx = 0
        prev = ''
        for c_idx, char in enumerate(wkt):
            if char == ',' and prev != ')':
                if insert_idx == point_idx:
                    insert_string = f', {pos[0]} {pos[1]}'
                    wkt = wkt[:c_idx] + insert_string + wkt[c_idx:]
                    break
                point_idx += 1
            prev = char
        _, bin = self.compress(shapely.from_wkt(wkt))