#!/usr/bin/env python3

'''
/**
 * @file escher_pattern.py
 * @brief wrapper for web
 * @par <b>License</b>:
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
'''

from urllib.parse import urlparse, parse_qs

import time

# this is needed for importing?
import os, sys
os.chdir(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__))
# also add:
'''
/etc/apache2/conf-enabled/serve-cgi-bin.conf:
                <Directory /var/www/html/pytest>
                        Options +ExecCGI
                        AddHandler wsgi-script .py
                        WSGIProcessGroup %{GLOBAL}
                        WSGIApplicationGroup %{GLOBAL}
                        Order deny,allow
                        Allow from all
                </Directory>

'''

from escher_pattern_class import Escher_Pattern

# def application is a hardcoded keyword for WSGI ?!
def application(environ,start_response):

    params = parse_qs(environ["QUERY_STRING"])

    status = '200 OK'

    lpm = 50
    page_width  = 270
    page_height = 210
    rotate = 5
    escher = 2.0
    tmp_dir = 'tmp'

    try:
      lpm = float(params['LPM'][0])
    except KeyError:
      pass

    try:
      page_width = float(params['PAGE_WIDTH'][0])
    except KeyError:
      pass

    try:
      page_height = float(params['PAGE_HEIGHT'][0])
    except KeyError:
      pass

    try:
      rotate = float(params['ROTATE'][0])
    except KeyError:
      pass

    try:
      escher = float(params['ESCHER'][0])
    except KeyError:
      pass

    basename  = 'escher-pattern'
    basename += '-ESCHER'+str(escher)
    basename += '-LPM'+str(lpm)
    basename += '-ROT'+str(rotate)
    basename += '-PAGE_WIDTH'+str(page_width)
    basename += '-PAGE_HEIGHT'+str(page_height)
    pdf_name  = basename+".pdf"

    ep = Escher_Pattern(width= page_width, height= page_height, escher=escher, lpm=lpm, rotate=rotate)
    ep.generate()
    ep.save()

    with open(pdf_name,'rb') as f:
      contents = f.read()

    filesize = str(len(contents))
    filesize = filesize.encode('utf-8')

    response_header = [
      ('Content-Type','application/pdf'),
      ('Cache-Control', 'public, must-revalidate, max-age=0'),
      ('Pragma', 'public'),
      ('Expires', 'Sat, 26 Jul 1997 05:00:00 GMT'),
      ('Last-Modified', time.strftime("%a, %d %b %Y %H:%M:%S GMT", time.gmtime())),
      ('Content-Length', str(len(contents))),
      ('Content-Disposition', 'inline; filename='+pdf_name+';')
    ]

    start_response(status,response_header)

    return [contents]
