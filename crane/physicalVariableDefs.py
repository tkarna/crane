"""
Definitions of variable names, units etc.

Tuomas Karna 2015-08-21
"""

# Long human readable names for variables
# Used in figure labels etc.

VARS = {'salt': 'Salinity',
        'temp': 'Temperature',
        'cond': 'Electrical conductivity',
        'elev': 'Elevation',
        'pres': 'Pressure',
        'fluores': 'Fluorescence',
        'kine': 'TKE',
        'vdff': 'Viscosity',
        'tdff': 'Diffusivity',
        'mixl': 'Mixing length',
        'hvel': 'Velocity',
        'vert': 'Vert velocity',
        'u': 'u-velocity',
        'v': 'v-velocity',
        'w': 'w-velocity',
        'alongvel': 'Along Vel.',
        'discharge': 'Discharge',
        'tidal_range': 'Tidal range',
        'cui': 'Upwelling index',
        'strat': 'Stratification',
        'sil': 'Salt intrusion length',
        'srad': 'Solar radiation',
        'sed': 'Sediment',
        'bed_stress': 'bed stress',
        'bed_depth': 'bed depth',
        'sed_1': 'S-class 1',
        'sed_2': 'S-class 2',
        'sed_3': 'S-class 3',
        'sed_4': 'S-class 4',
        'di_sed': 'DI sed',
        'NO3': 'NO3',
        'NH4': 'NH4',
        'PHYM': 'PHYM',
        'PHYF': 'PHYF',
        'SZOO': 'SZOO',
        'BZOO': 'BZOO',
        'DETN': 'DETN',
        'DETC': 'DETC',
        'BACT': 'BACT',
        'DON': 'DON',
        'DOC': 'DOC',
        'CHLM': 'CHLM',
        'CHLF': 'CHLF',
        'bnth_1': 'NO3 benthos',
        'bnth_2': 'NH4 benthos',
        'bnth_3': 'DETN benthos',
        'bnth_4': 'OXY benthos',
        'Diag': 'Diagnostic array',
        'nem': 'NEM',
        'nem_DAVG': 'Averaged NEM',
        'nem_DI': 'Integrated NEM',
        'nemb': 'benthic NEM',
        'nemi': 'vert. int. NEM',
        'totalN': 'total Nitrogen',
        'NO3': 'NO3',
        'NH4': 'NH4',
        'phy': 'Phy',
        'zoo': 'Zoo',
        'det': 'Det',
        'oxy': 'Oxy',
        'plumevol': 'Plume volume',
        'plumethick': 'Plume thickness',
        'flntu': 'FLNTU',
        'ecoblue': 'ECO backscatter 470nm',
        'ecored': 'ECO backscatter 700nm',
        'dens': 'Density',
        'turbidity': 'Sediment',
        'bath': 'Bathymetry',
        'depth': 'Depth',
        'bed_depth': 'Depth',
        'bed_stress': 'bottom stress',
        'sil_1': 'SIL 1',
        'sil_2': 'SIL 2',
        'sil_3': 'SIL 3',
        'sil_4': 'SIL 4',
        'sil_5': 'SIL 5',
        'plume_area_12': 'Plume area 12',
        'plume_volume_12': 'Plume volume 12',
        'plume_center_12': 'Plume center 12',
        'plume_center_x_12': 'Plume center x 12',
        'plume_center_y_12': 'Plume center y 12',
        'plume_thickness_12': 'Plume thickness 12',
        'plume_area_28': 'Plume area 28',
        'plume_volume_28': 'Plume volume 28',
        'plume_center_28': 'Plume center 28',
        'plume_center_x_28': 'Plume center x 28',
        'plume_center_y_28': 'Plume center y 28',
        'plume_thickness_28': 'Plume thickness 28',
        'plume_area_30': 'Plume area 30',
        'plume_volume_30': 'Plume volume 30',
        'plume_center_30': 'Plume center 30',
        'plume_center_x_30': 'Plume center x 30',
        'plume_center_y_30': 'Plume center y 30',
        'plume_thickness_30': 'Plume thickness 30',
        'volume': 'Volume',
        'flux': 'Flux',
        'bed_flux': 'Sediment bed flux',
        'iRiv': 'Riverine water',
        'iOce': 'Oceanic water',
        'iPlu': 'Plume water',
        'iRen': 'Renewed water',
        'aRiv': 'Riv. age conc.',
        'aOce': 'Oce. age conc.',
        'aPlu': 'Plu. age conc.',
        'aRen': 'Ren. age conc.',
        'ARiv': 'Riv. age',
        'AOce': 'Oce. age',
        'APlu': 'Plu. age',
        'ARen': 'Ren. age',
        'maxiRiv': 'Riverine water',
        'maxiOce': 'Oceanic water',
        'maxiPlu': 'Plume water',
        'maxiRen': 'Renewed water',
        'maxaRiv': 'Riv. age conc.',
        'maxaOce': 'Oce. age conc.',
        'maxaPlu': 'Plu. age conc.',
        'maxared': 'Ren. age conc.',
        'maxARiv': 'Riv. age',
        'maxAOce': 'Oce. age',
        'maxAPlu': 'Plu. age',
        'maxARen': 'Ren. age',
        'daviRiv': 'Depth av. Riverine water',
        'daviOce': 'Depth av. Oceanic water',
        'daviPlu': 'Depth av. Plume water',
        'daviRen': 'Depth av. Renewed water',
        'davaRiv': 'Depth av. Riv. age conc.',
        'davaOce': 'Depth av. Oce. age conc.',
        'davaPlu': 'Depth av. Plu. age conc.',
        'davared': 'Depth av. Ren. age conc.',
        'davARiv': 'Depth av. Riverine water age',
        'davAOce': 'Depth av. Oceanic water age',
        'davAPlu': 'Depth av. Plume water age',
        'davARen': 'Depth av. renewal water age',
        'tavdaviRiv': 'Riverine water',
        'tavdaviOce': 'Oceanic water',
        'tavdaviPlu': 'Plume water',
        'tavdaviRen': 'Renewed water',
        'tavdavaRiv': 'Riv. age conc.',
        'tavdavaOce': 'Oce. age conc.',
        'tavdavaPlu': 'Plu. age conc.',
        'tavdavared': 'Ren. age conc.',
        'tavdavARiv': 'Riverine water age',
        'tavdavAOce': 'Oceanic water age',
        'tavdavAPlu': 'Plume water age',
        'tavdavARen': 'renewal water age',
        'vort': 'Vorticity',
        'tavvort': 'Time av. Vorticity',
        'davvort': 'Depth av. Vorticity',
        'tavdavvort': 'Vorticity',
        }

# Units for variables
# Used in figure labels etc.
UNITS = {'salt': 'psu',
         'temp': 'degC',
         'cond': '',
         'elev': 'm',
         'pres': 'Pa',
         'fluores': '',
         'kine': 'm^2/s^2',
         'vdff': 'm^2/s',
         'tdff': 'm^2/s',
         'mixl': 'm',
         'hvel': 'm/s',
         'vert': 'm/s',
         'u': 'm/s',
         'v': 'm/s',
         'w': 'm/s',
         'alongvel': 'm/s',
         'discharge': 'm^3/s',
         'tidal_range': 'm',
         'cui': '',
         'strat': 'psu',
         'sil': 'km',
         'srad': 'W m-2',
         'sed': 'mg/l',
         'bed_stress': 'N/m^2',
         'bed_depth': 'm',
         'sed_1': 'kg/m^3',
         'sed_2': 'kg/m^3',
         'sed_3': 'kg/m^3',
         'sed_4': 'kg/m^3',
         'sed_5': 'kg/m^3',
         'di_sed': 'kg/m^3',
         'NO3': 'mmol N m-3',
         'NH4': 'mmol N m-3',
         'PHYM': 'mmol N m-3',
         'PHYF': 'mmol N m-3',
         'SZOO': 'mmol N m-3',
         'BZOO': 'mmol N m-3',
         'DETN': 'mmol N m-3',
         'DETC': 'mmol C m-3',
         'BACT': 'mmol N m-3',
         'DON': 'mmol N m-3',
         'DOC': 'mmol C m-3',
         'CHLM': 'mg Chl m-3',
         'CHLF': 'mg Chl m-3',
         'bnth_1': 'mmol N m-3',
         'bnth_2': 'mmol N m-3',
         'bnth_3': 'mmol N m-3',
         'bnth_4': 'mmol O m-3',
         'Diag': ' ',
         'nem': 'mmol O m-3 day-1 ',
         'nem_DAVG': 'mmol O m-3 day-1 ',
         'nem_DI': 'mmol O m-2 day-1 ',
         'nemb': 'mmol O m-2 day-1 ',
         'nemi': 'mmol O m-2 day-1 ',
         'totalN': 'mmol N m-3',
         'trcr_5': 'kg/m^3',
         'phy': 'mmol N m-3',
         'zoo': 'mmol N m-3',
         'det': 'mmol N m-3',
         'oxy': 'mmol m-3',
         'plumevol': '10^9 m^3',
         'plumethick': 'm',
         'flntu': 'ntu',
         'ecored': 'm-1 sr-1',
         'ecoblue': 'm-1 sr-1',
         'dens': 'kg m-3',
         'turbidity': 'kg m-3',
         'SPM': 'kg m-3',
         'bath': 'm',
         'depth': 'm',
         'sil_1': 'km',
         'sil_2': 'km',
         'sil_3': 'km',
         'sil_4': 'km',
         'sil_5': 'km',
         'plume_area_12': 'm2',
         'plume_volume_12': 'm3',
         'plume_center_12': 'm',
         'plume_center_x_12': 'm',
         'plume_center_y_12': 'm',
         'plume_thickness_12': 'm',
         'plume_area_28': 'm2',
         'plume_volume_28': 'm3',
         'plume_center_28': 'm',
         'plume_center_x_28': 'm',
         'plume_center_y_28': 'm',
         'plume_thickness_28': 'm',
         'plume_area_30': 'm2',
         'plume_volume_30': 'm3',
         'plume_center_30': 'm',
         'plume_center_x_30': 'm',
         'plume_center_y_30': 'm',
         'plume_thickness_30': 'm',
         'volume': 'm3',
         'flux': 'm3 s-1',
         'bed_flux': 'm3',
         'iRiv': 'frac',
         'iOce': 'frac',
         'iPlu': 'frac',
         'iRen': 'frac',
         'aRiv': 'frac s',
         'aOce': 'frac s',
         'aPlu': 'frac s',
         'aRen': 'frac s',
         'ARiv': 'h',
         'AOce': 'h',
         'APlu': 'h',
         'ARen': 'h',
         'maxiRiv': 'frac',
         'maxiOce': 'frac',
         'maxiPlu': 'frac',
         'maxiRen': 'frac',
         'maxaRiv': 'frac s',
         'maxaOce': 'frac s',
         'maxaPlu': 'frac s',
         'maxaRen': 'frac s',
         'maxARiv': 'h',
         'maxAOce': 'h',
         'maxAPlu': 'h',
         'maxARen': 'h',
         'daviRiv': 'frac',
         'daviOce': 'frac',
         'daviPlu': 'frac',
         'daviRen': 'frac',
         'davaRiv': 'frac s',
         'davaOce': 'frac s',
         'davaPlu': 'frac s',
         'davaRen': 'frac s',
         'davARiv': 'h',
         'davAOce': 'h',
         'davAPlu': 'h',
         'davARen': 'h',
         'tavdaviRiv': 'frac',
         'tavdaviOce': 'frac',
         'tavdaviPlu': 'frac',
         'tavdaviRen': 'frac',
         'tavdavaRiv': 'frac s',
         'tavdavaOce': 'frac s',
         'tavdavaPlu': 'frac s',
         'tavdavaRen': 'frac s',
         'tavdavARiv': 'h',
         'tavdavAOce': 'h',
         'tavdavAPlu': 'h',
         'tavdavARen': 'h',
         'vort': 's-1',
         'tavvort': 's-1',
         'davvort': 's-1',
         'tavdavvort': 's-1',
         }

# Associates each variable to a SELFE output file
fieldNameToFilename = {'temp': 'temp.63',
                       'elev': 'elev.61',
                       'salt': 'salt.63',
                       'alongvel': 'hvel.64',
                       'kine': 'kine.63',
                       'vdff': 'vdff.63',
                       'tdff': 'tdff.63',
                       'mixl': 'mixl.63',
                       'hvel': 'hvel.64',
                       'vert': 'vert.63',
                       'dens': 'conc.63',
                       'srad': 'srad.61',
                       'dahv': 'dahv.62',
                       'wind': 'wind.62',
                       'wist': 'wist.62',
                       'vort': 'vort.63',
                       'trcr_1': 'trcr_1.63',
                       'trcr_2': 'trcr_2.63',
                       'trcr_3': 'trcr_3.63',
                       'trcr_4': 'trcr_4.63',
                       'trcr_5': 'trcr_5.63',
                       'trcr_6': 'trcr_6.63',
                       'trcr_7': 'trcr_7.63',
                       'trcr_8': 'trcr_8.63',
                       'trcr_9': 'trcr_9.63',
                       'trcr_10': 'trcr_10.63',
                       'trcr_11': 'trcr_11.63',
                       'trcr_12': 'trcr_12.63',
                       'trcr_13': 'trcr_13.63',
                       'trcr_14': 'trcr_14.63',
                       'bnth_1': 'bnth_1.61',
                       'bnth_2': 'bnth_2.61',
                       'bnth_3': 'bnth_3.61',
                       'bnth_4': 'bnth_4.61',
                       'carai': 'nem_DAVG.61',
                       'Diag': 'Diag.63',
                       'Dia2': 'Dia2.61',
                       'nem': 'nem.63',
                       'nem_DAVG': 'nem_DAVG.61',
                       'nem_DI': 'nem_DI.61',
                       'nemb': 'nemb.61',
                       'nemi': 'nemi.61',
                       'prod': 'prod.63',
                       'prfr': 'prfr.63',
                       'prma': 'prma.63',
                       'resp': 'resp.63',
                       'rezo': 'rezo.63',
                       'reba': 'reba.63',
                       'ream': 'ream.63',
                       'rere': 'rere.63',
                       'totalN': 'totalN.63',
                       'turbidity': 'turbidity.63',
                       'bed_age': 'bed_age.61',
                       'bed_depth': 'bed_depth.61',
                       'bed_mass': 'bed_mass.61',
                       'bed_thick': 'bed_thick.61',
                       'bed_stress': 'bed_stress.62',
                       'bed_load_1': 'bed_load_1.62',
                       'bed_load_2': 'bed_load_2.62',
                       'bed_load_3': 'bed_load_3.62',
                       'bed_load_4': 'bed_load_4.62',
                       'bed_flux': 'bed_flux.61'}

# List of components of variables (e.g. velocity)

fieldNameList = {'temp': ['temp'],
                 'elev': ['elev'],
                 'salt': ['salt'],
                 'alongvel': ['u', 'v'],
                 'kine': ['kine'],
                 'vdff': ['vdff'],
                 'tdff': ['tdff'],
                 'mixl': ['mixl'],
                 'hvel': ['u', 'v'],
                 'wind': ['u', 'v'],
                 'wist': ['u', 'v'],
                 'vert': ['w'],
                 'dens': ['dens'],
                 'dahv': ['u', 'v'],
                 'trcr_1': ['trcr_1'],
                 'trcr_2': ['trcr_2'],
                 'trcr_3': ['trcr_3'],
                 'trcr_4': ['trcr_4'],
                 'trcr_5': ['trcr_5'],
                 'trcr_6': ['trcr_6'],
                 'trcr_7': ['trcr_7'],
                 'trcr_8': ['trcr_8'],
                 'trcr_9': ['trcr_9'],
                 'trcr_10': ['trcr_10'],
                 'trcr_11': ['trcr_11'],
                 'trcr_12': ['trcr_12'],
                 'trcr_13': ['trcr_13'],
                 'trcr_14': ['trcr_14'],
                 'bnth_1': ['bnth_1'],
                 'bnth_2': ['bnth_2'],
                 'bnth_3': ['bnth_3'],
                 'bnth_4': ['bnth_4'],
                 'carai': ['nem_DAVG'],
                 'Diag': ['Diag'],
                 'Dia2': ['Dia2'],
                 'nem': ['nem'],
                 'nem_DAVG': ['nem_DAVG'],
                 'nem_DI': ['nem_DI'],
                 'nemb': ['nemb'],
                 'nemi': ['nemi'],
                 'prod': ['prod'],
                 'prfr': ['prfr'],
                 'prma': ['prma'],
                 'resp': ['resp'],
                 'rezo': ['rezo'],
                 'reba': ['reba'],
                 'ream': ['ream'],
                 'rere': ['rere'],
                 'totalN': ['totalN'],
                 'turbidity': ['turbidity'],
                 'bed_age': ['bed_age'],
                 'bed_depth': ['bed_depth'],
                 'bed_mass': ['bed_mass'],
                 'bed_thick': ['bed_thick'],
                 'bed_stress': ['bed_stress'],
                 'bed_load_1': ['bed_load_1'],
                 'bed_load_2': ['bed_load_2'],
                 'bed_load_3': ['bed_load_3'],
                 'bed_load_4': ['bed_load_4'],
                 'bed_flux': ['bed_flux']}

# observational data that can be compared to model outputs for each tracer
# model
tracerModelObsVariables = {'sed': ['turbidity'],
                           'oxy': ['NO3', 'oxy'],
                           'bio': ['NO3', 'oxy']
                           }

#-------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------


def addTracers(tracerModelName, varList=None, numTracers=None):
    """Appends constants used throughout processing to support tracers from coupled
    models.
    """
    print 'Adding tracer model ' + tracerModelName
    global fieldNameToFilename
    global fieldNameList

    if varList is None:
        varList = []
    else:
        varList = list(varList)

    # parse tracerModelName 'oxy' or 'oxy.70'
    if len(tracerModelName.split('.')) == 2:
        tracerName, fileExtension = tracerModelName.split('.')
    else:
        tracerName = tracerModelName
        fileExtension = '63'

    if tracerName == 'oxy':
        fieldNameToFilename['NO3'] = 'trcr_1.' + fileExtension
        fieldNameToFilename['NH4'] = 'trcr_2.' + fileExtension
        fieldNameToFilename['phy'] = 'trcr_3.' + fileExtension
        fieldNameToFilename['zoo'] = 'trcr_4.' + fileExtension
        fieldNameToFilename['det'] = 'trcr_5.' + fileExtension
        fieldNameToFilename['oxy'] = 'trcr_6.' + fileExtension

        fieldNameList['NO3'] = ['NO3']
        fieldNameList['NH4'] = ['NH4']
        fieldNameList['phy'] = ['Phy']
        fieldNameList['zoo'] = ['Zoo']
        fieldNameList['det'] = ['Det']
        fieldNameList['oxy'] = ['Oxy']

        trcr_vars = ['NO3', 'NH4', 'phy', 'zoo', 'det', 'oxy']

        for i in range(1, len(trcr_vars) + 1):
            varList.append('trcr_%d.%s' % (i, fileExtension))

    elif tracerName == 'oxy2':
        fieldNameToFilename['oxy'] = 'trcr_1.' + fileExtension
        fieldNameList['oxy'] = ['oxygen']
        trcr_vars = ['oxy']
        varList.append('trcr_1.%s' % (fileExtension))

    elif tracerName == 'bio':
        # these are internal names used in processing lib to identify tracers
        # these can be anything, they appear in file names, but not in plots
        trcr_vars = ['NO3', 'NH4', 'PHYM', 'PHYF', 'SZOO', 'BZOO', 'DETN',
                     'DETC', 'BACT', 'DON', 'DOC', 'CHLM', 'CHLF', 'oxy']
        # maps each var to selfe output file
        for i, v in enumerate(trcr_vars):
            fieldNameToFilename[v] = 'trcr_{0:d}.{1:s}'.format(
                i + 1, fileExtension)

        # these will be stored in dataContainer fieldnames list
        # to distinguish different components (e.g. hvel -> [u,v])
        # for tracers the same string can be used
        for v in trcr_vars:
            fieldNameList[v] = [v]

        for i in range(len(trcr_vars)):
            varList.append('trcr_{0:d}.{1:s}'.format(i + 1, fileExtension))

    elif tracerName == 'age':
        # water age model with 2 indicator tracers: river and ocean
        trcr_vars = ['iRiv', 'iOce', 'aRiv', 'aOce']
        for i, v in enumerate(trcr_vars):
            fieldNameToFilename[v] = 'trcr_{0:d}.{1:s}'.format(
                i + 1, fileExtension)

        for v in trcr_vars:
            fieldNameList[v] = [v]

        for i in range(len(trcr_vars)):
            varList.append('trcr_{0:d}.{1:s}'.format(i + 1, fileExtension))

    elif tracerName == 'age2':
        # water age model with 3 indicator tracers: river, ocean and plume
        trcr_vars = ['iRiv', 'iOce', 'iPlu', 'aRiv', 'aOce', 'aPlu',
                     'ARiv', 'AOce', 'APlu']
        for i, v in enumerate(trcr_vars):
            fieldNameToFilename[v] = 'trcr_{0:d}.{1:s}'.format(
                i + 1, fileExtension)

        for v in trcr_vars:
            fieldNameList[v] = [v]

        for i in range(len(trcr_vars)):
            varList.append('trcr_{0:d}.{1:s}'.format(i + 1, fileExtension))

    elif tracerName == 'sed':
        for t in range(1, numTracers + 1):
            fieldNameToFilename['sed_%d' %
                                t] = 'trcr_%d.%s' % (t, fileExtension)
            fieldNameList['sed_%d' % t] = ['sed_class_%d' % t]
            varList.append('trcr_%d.%s' % (t, fileExtension))
    elif tracerName == 'generic':
        for t in range(1, numTracers + 1):
            fieldNameToFilename['trcr_%d' %
                                t] = 'trcr_%d.%s' % (t, fileExtension)
            fieldNameList['trcr_%d' % t] = ['trcr_%d' % t]
            varList.append('trcr_%d.%s' % (t, fileExtension))
    else:
        raise Exception('Tracer model not supported: ' + tracerName)

    return varList
