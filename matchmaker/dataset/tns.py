from ...extern.tns_api_search.tns_api_search import search, format_to_json
import urllib.request
import pandas as pd
import requests
from bs4 import BeautifulSoup
from astropy.coordinates import SkyCoord
import time


# This should go somewhere safe.
TNS_CREDENTIALS = {
    'TNS_BOT_ID': '131254',
    'TNS_BOT_NAME': 'matchmaker',
    'TNS_API_KEY': '4297df75861837b6260c1f4734ef0643269bae81'
}
'''
row2json, get_json, json2df 
     taken/modified from 
     https://github.com/davidgardenier/frbcat/tns.py
'''
def row2json(line):
    """Convert row of html table to json format."""
    out = {}
    for seg in line.split('</td>'):
        if '<td ' in seg:
            key = seg.split('<td ')[-1].split('class="cell-')[-1]
            key = key.split('"')[0]
            val = seg.split('<td ')[-1].split('>', 1)[1].strip()
            if val:
                ks = ('filename', 'public_webpage', 'region_filename')
                if key in ks:
                    val = val.split('href="')[1].split('"')[0]
                if key in ('photometry', 'related_files', 'reps'):
                    val = val.split('<a')[0]
                if key in ('id', 'name', 'repeater_of_objid'):
                    val = val.split('</a>')[0].split('>')[-1]
                out[key] = val
            # if '<' in val:
            #     print(key, val)
    return out

def get_results_as_dataframe(page = None):
    more = True
    page = 1 if page is None else page
    pause = 0
    num_page = 500
    data = []
    columns = []

    # Provide user agent to be able to access the webpage
    header = {'User-Agent': str({
        'tns_id': TNS_CREDENTIALS['TNS_BOT_ID'],
         'type': 'user',
         'name': TNS_CREDENTIALS['TNS_BOT_NAME']
    })}

    # Loop through pages on TNS webpage till no more results
    print ('Start fetching TNS pages')
    while more:
        # Limit results classified SN
        url = "https://www.wis-tns.org/search?" \
        "&discovered_period_value=" \
        "&discovered_period_units=days" \
        "&unclassified_at=0" \
        "&classified_sne=1" \
        "&include_frb=0" \
        "&name=" \
        "&name_like=0" \
        "&isTNS_AT=all" \
        "&public=all" \
        "&ra=" \
        "&decl=" \
        "&radius=" \
        "&coords_unit=arcsec" \
        "&reporting_groupid%5B%5D=null" \
        "&groupid%5B%5D=null" \
        "&classifier_groupid%5B%5D=null" \
        "&objtype%5B%5D=null" \
        "&at_type%5B%5D=null" \
        "&date_start%5Bdate%5D=" \
        "&date_end%5Bdate%5D=" \
        "&discovery_mag_min=" \
        "&discovery_mag_max=" \
        "&internal_name=" \
        "&discoverer=" \
        "&classifier=" \
        "&spectra_count=" \
        "&redshift_min=" \
        "&redshift_max=" \
        "&hostname=" \
        "&ext_catid=" \
        "&ra_range_min=" \
        "&ra_range_max=" \
        "&decl_range_min=" \
        "&decl_range_max=" \
        "&discovery_instrument%5B%5D=null" \
        "&classification_instrument%5B%5D=null" \
        "&associated_groups%5B%5D=null" \
        "&official_discovery=0" \
        "&official_classification=0" \
        "&at_rep_remarks=" \
        "&class_rep_remarks=" \
        "&frb_repeat=all" \
        "&frb_repeater_of_objid=" \
        "&frb_measured_redshift=0" \
        "&frb_dm_range_min=" \
        "&frb_dm_range_max=" \
        "&frb_rm_range_min=" \
        "&frb_rm_range_max=" \
        "&frb_snr_range_min=" \
        "&frb_snr_range_max=" \
        "&frb_flux_range_min=" \
        "&frb_flux_range_max=" \
        "&num_page={}" \
        "&display%5Bredshift%5D=1" \
        "&display%5Bhostname%5D=1" \
        "&display%5Bhost_redshift%5D=1" \
        "&display%5Bsource_group_name%5D=1" \
        "&display%5Bclassifying_source_group_name%5D=1" \
        "&display%5Bdiscovering_instrument_name%5D=0" \
        "&display%5Bclassifing_instrument_name%5D=0" \
        "&display%5Bprograms_name%5D=0" \
        "&display%5Binternal_name%5D=1" \
        "&display%5BisTNS_AT%5D=0" \
        "&display%5Bpublic%5D=1" \
        "&display%5Bend_pop_period%5D=0" \
        "&display%5Bspectra_count%5D=1" \
        "&display%5Bdiscoverymag%5D=1" \
        "&display%5Bdiscmagfilter%5D=1" \
        "&display%5Bdiscoverydate%5D=1" \
        "&display%5Bdiscoverer%5D=1" \
        "&display%5Bremarks%5D=0" \
        "&display%5Bsources%5D=0" \
        "&display%5Bbibcode%5D=0" \
        "&display%5Bext_catalogs%5D=0" \
        "{}".format(num_page, '&page={}'.format(page) if page > 1 else '')

        # print (url)
        # url = "https://www.wis-tns.org/search?" \
        #       "&page={}" \
        #       "&unclassified_at=0" \
        #       "&classified_sne=1" \
        #       "&include_frb=0" \
        #       "&name_like=0" \
        #       "&isTNS_AT=all" \
        #       "&public=all" \
        #       "&num_page={}" \
        #       "&format=html" \
        #       "&edit[type]=" \
        #       "&edit[objname]=" \
        #       "&edit[id]=".format(page, num_page)

        t = time.process_time()
        html = requests.get(url, headers=header)
        soup = BeautifulSoup(html.text)
        soup.find('li', {'class':'pager-item'})
        table = soup.find('table', {'class':'results-table'})
        if len(columns) == 0:
            columns = [i.text for i in table.thead.find_all('th')]
            n_pages = int(str(soup.find('li', {'class':'pager-last'}).find('a')).split('?page=')[-1].split('&')[0])
        # Get data
        tmp = [[cell.text for cell in row.find_all('td')] for row in table.tbody.select('tr.public.odd')]
        for d in tmp:
            # rows on the website include hidden rows. The main rows have 29 items specifically. This is a bit of a hack
            if len(d) == 29:
                data.append(d)

        # print (page, len(data), (num_page * (page+1)))
        # if len(data) > 0 and len(data) % (num_page * (page+1)) == 0:


        print ('Fetched %d out of %d pages. %.2f sec' % (page, n_pages, time.process_time() - t))
        if page < n_pages:
            page += 1
            # if page % 14 == 0:
            #     print ('try wait for 60 seconds to avoid being bumped')
            #     time.sleep(60)
        else:
            more = False
            print ('Fetch completed')


    df = pd.DataFrame(data, columns=columns)
    df['ra'] = df.apply(lambda row: SkyCoord('{}h{}m{}s'.format(row['RA'].split(':')[0], row['RA'].split(':')[1], row['RA'].split(':')[2]),
                                           '{}d{}m{}s'.format(row['DEC'].split(':')[0], row['DEC'].split(':')[1], row['DEC'].split(':')[2]), frame='icrs').ra.deg, axis=1)
    df['dec'] = df.apply(lambda row: SkyCoord('{}h{}m{}s'.format(row['RA'].split(':')[0], row['RA'].split(':')[1], row['RA'].split(':')[2]),
                                                '{}d{}m{}s'.format(row['DEC'].split(':')[0], row['DEC'].split(':')[1], row['DEC'].split(':')[2]), frame='icrs').dec.deg, axis=1)

    # Coerce removes anything before pd.Timestamp.min
        # namely, 1677-09-21 00:12:43.145224193
        # Otherwise, pandas raises an error
    df['Discovery Date (UT)'] = pd.to_datetime(df['Discovery Date (UT)'], errors = 'coerce', utc=True)
    df['Class'] = pd.to_numeric(df['Class'])
    df['Redshift'] = pd.to_numeric(df['Redshift'])
    df['Discovery Mag/Flux'] = pd.to_numeric(df['Discovery Mag/Flux'])

    return df

def json2df(entries):
#
    # Create a nice list of dictionaries
    rows = []
    for frb in entries:
        row = {}
        for par in frb:
            if type(frb[par]) == list and len(frb[par]) > 0:
                # Always take the most recent entry
                # TODO: This might need to be updated at some stage
                for other_par in frb[par][0]:
                    name = other_par
                    if other_par in frb:
                        name = par.split('_')[0] + '_' + other_par
                    row[name] = frb[par][0][other_par]
            else:
                row[par] = frb[par]
        rows.append(row)

    # Convert to a DataFrame
    return pd.DataFrame(rows)

def query():
    search_obj = [
        ("classified_sne", 1),
        ("objname", ""),
        ("objname_exact_match", 0),
        ("internal_name", ""),
        ("internal_name_exact_match", 0),
        ("objid", ""),
        ("public_timestamp", "")
    ]
    response = search(search_obj, credentials=TNS_CREDENTIALS)
    json_data = format_to_json(response.text)



    return response, json_data, search_obj
