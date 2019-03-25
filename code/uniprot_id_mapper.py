def uniprot_id_mapper(ids, origindb, targetdb ):
    import urllib
    import math
    import pandas as pd
    from pandas.compat import StringIO

    url = 'https://www.uniprot.org/uploadlists/'
    id_mapping=pd.DataFrame()

    # 10 slices
    slices = (ids[(i*10):min(len(ids),(i+1)*10)] for i in range(math.ceil(len(ids)/10)))

    for id_slice in slices:
        params = {
        'from':origindb,
        'to':targetdb,
        'format':'tab',
        'query':id_slice
        }

        data = urllib.parse.urlencode(params)
        data = data.encode('ascii')
        request = urllib.request.Request(url, data=data)
        contact = "simon.couvreur@kcl.ac.uk" # Please set a contact email address here to help us
            # debug in case of problems (see https://www.uniprot.org/help/privacy).
        request.add_header('User-Agent', 'Python %s' % contact)
        response = urllib.request.urlopen(request)
        slice_map = pd.read_csv(StringIO(response.read(200000).decode()), sep="\t")
        # print(slice_map)
        if not slice_map.empty:
            id_mapping = id_mapping.append(slice_map)
    return(id_mapping)
