#!/usr/bin/env python
import json
import pycurl, cStringIO
buf = cStringIO.StringIO()
data = {
    'token': 'B7135B24799ECCCF7B80903AA3D7202A',
    'content': 'project',
    'format': 'json',
    'returnFormat': 'json'
}
ch = pycurl.Curl()
ch.setopt(ch.URL, 'http://10.67.1.22/redcap/api/')
ch.setopt(ch.HTTPPOST, data.items())
ch.setopt(ch.WRITEFUNCTION, buf.write)
ch.perform()
ch.close()
jsondata = json.loads(buf.getvalue())
buf.close()

json_text = json.dumps(jsondata,indent=4)
print json_text

# list of content
# instrument : liste des instruments (sans details)
# exportFieldNames : liste des champs des instruments, nom et choix possibles
# metadata : metadonnees sur la maniere dont sont construits les differents champs
