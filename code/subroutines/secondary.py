# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 12:46:24 2020

@author: aargl
"""

import os
import numpy as np
import pandas as pd
import datetime

#%% Parsing Data
    

def load_COVID19_mort():
    
    """
    UK release of daily deaths arrising from likely COVID-19 cases in hosptials
    and care homes.
    
    (Accessed: 19/06/2020) https://coronavirus.data.gov.uk/
    
    Inputs:
        
        None
        
    Outputs
    
        - date, the date of individual data points.
        
        - mort, the mortality rate. (units: deaths)
        
        - R, the reproduction number. (units: )
    
    """
    

    # Load in Data
    filename = 'coronavirus-deaths_latest.xlsx'
    excel_spreadsheet =  loaddata(filename,domain='covid19')
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[1]) 
    
    # Get Date
    date_str = np.array(df.iloc[0:106,3].astype(str))
    num_days = len(date_str)
    date = np.zeros(num_days,dtype=datetime.datetime)
    
    # Load date
    for day in range(0,num_days):
        
        date[day] = datetime.datetime.strptime(date_str[day].split('T')[0],'%Y-%m-%d')
        
    # Load in daily death rate
    mort = np.array(df.iloc[0:106,4])
    # Estimate the reproduction rate 
    
    R = np.nan*np.ones(num_days)
    for day in range(1,num_days):
        
        if mort[day-1] > 0:
            
            R[day] = mort[day]/mort[day-1]

    return date, mort, R

def load_google_mobility():
    
    """
    Loads in the mobility data from Google for the United Kingdom. Interprets a
    general m from the daily maximum of either ‘transit stations’, ‘workplaces’ and ‘retail and recreation’ mobilities
    
    Google LLC "Google COVID-19 Community Mobility Reports".
    https://www.google.com/covid19/mobility/ Accessed: 18/05/2020
    
    Inputs:
        
        None
        
    Outputs
    
        - date, the date of individual data points 
        
        - m, the mobility (units: fraction change from baseline)
    
    """
    

    # Load in Data
    filename = 'UK_Mobility.xlsx'
    excel_spreadsheet =  loaddata(filename,domain='covid19')
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[0]) 
    
    # Get Date
    date_str = np.array(df.iloc[0:120,6].astype(str))
    num_days = len(date_str)
    date = np.zeros(num_days,dtype=datetime.datetime)
    
    # Load date
    for day in range(0,num_days):
        
        date[day] = datetime.datetime.strptime(date_str[day].split('T')[0],'%Y-%m-%d')
        
    # Load "retail and recreataion", "workplaces" and "transit stations" rate of change of mobility from baseline
    r_m_rr = np.array(df.iloc[0:120,7])
    r_m_wp = np.array(df.iloc[0:120,10])
    r_m_ts = np.array(df.iloc[0:120,11])
    # convert to relative to baseline
    m_3 = np.zeros((3,num_days))
    m_3[0,:] =  1+(r_m_rr)/100
    m_3[1,:] = 1+(r_m_wp)/100
    m_3[2,:] = 1+(r_m_ts)/100    
    m = np.nanmax(m_3,axis=0)
    
    return date, m
    

def load_covid_GDP_ons():
    
    """
    Load in Gross Domestic Product month on month change (GDP)
    
    Month Index from Jan 1997 using Gross Value Added from the Office for National Statistics: https://www.ons.gov.uk/economy/grossdomesticproductgdp/bulletins/gdpmonthlyestimateuk/april2020.
    
    Historical period is consistent with the ONS historical measurements.
    
    Inputs:
        
        None

    Outputs:
    
        - date, the year midpoint of the corresponding annual estimate. (type: datetime)
                
        - GDP, easonally adjusted Gross Domestic Product (GDP) (£billiion)
        
        - num_months, number of months of timeseries.
        
    
    """
    
    # Load in Data
    filename = 'GDP_change_ONS.xls'
    excel_spreadsheet =  loaddata(filename,domain='economy')
    # sheet_name[0] for month on month change to GDP
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[0])

    # Get Date
    date_str = np.array(df.iloc[6:287,0])
    num_months = len(date_str)
    date = np.zeros(num_months,dtype=datetime.datetime)    

    for mth in range(0,num_months):
        
        date[mth] = datetime.datetime.strptime(date_str[mth],'%Y %b')
        
        if date[mth].month < 12:

            date[mth] =  date[mth] + (datetime.datetime(date[mth].year,date[mth].month+1,1) - date[mth])/2 - datetime.timedelta(days=1)
            
        else:
            
            date[mth] =  date[mth] + (datetime.datetime(2020,1,1) - datetime.datetime(2019,12,1))/2 - datetime.timedelta(days=1)

    # Get fractional change of forecasted GDP
    rel_month = np.array(df.iloc[6:287,1]) # monthly index
    rel_month[0] = float(rel_month[0]) # first indicies was a string type for some reason
    
    r_GDP = rel_month[1:]/rel_month[0:num_months-1] - 1
    # Estimate in absolute terms, we assume end of the year in 1996 (£1305.527B) is equiv for jan + mean growth rate for that year https://www.ons.gov.uk/economy/grossdomesticproductgdp/timeseries/abmi/pn2
    gdp_1997 = 1355.853
    gdp_1996 = 1305.527
    gdp_hist_ons_mnth = np.zeros(num_months)
    gdp_hist_ons_mnth[0] = gdp_1996 + (gdp_1997/gdp_1996-1)/12
    
    for mnth in range(1,num_months):
        
        gdp_hist_ons_mnth[mnth] = gdp_hist_ons_mnth[mnth-1] * (1 + r_GDP[mnth-1])
        
    return date, gdp_hist_ons_mnth, num_months

def load_precovid_annual_GDP_estimate():
    
    """
    Load in government projection of seasonally adjusted Gross Domestic Product (GDP)
    
    Percent changes based on 2018 policies are given from the Department for Business, Energy & Industrial Strategy: https://www.gov.uk/government/publications/updated-energy-and-emissions-projections-2018
    
    Historical period is consistent with the ONS historical measurements.
    
    Inputs:
        
        None


    Outputs:
    
        - date, the year midpoint of the corresponding annual estimate. (type: datetime)
        
        - sim_elms, years that are simulated. (greater than 2016)
        
        - r_GDP, fractional change in GDP ()
        
        - num_years, number of years of time-series. 
        
    
    """
    
    # Load in Data
    filename = 'Annex-m-price-growth-assumption_16-May-2019.xlsx'
    excel_spreadsheet =  loaddata(filename,domain='economy')
    # sheet_name[1] for reference case used in emissions
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[1])
    # Get Year
    year = np.array(df.iloc[23,4:39])
    
    num_years = len(year)
    date = np.zeros(num_years,dtype=datetime.datetime)    

    for yr in range(0,num_years):
        
        date[yr] = datetime.datetime(int(year[yr]),1,1)  + (datetime.datetime(int(year[yr]),12,31) - datetime.datetime(int(year[yr]),1,1))/2.0
            
    num_years = len(year)
    
    # Year of simulation    
    sim_elms = date>datetime.datetime(2017,1,1)
    
    # Get fractional change of forecasted GDP (convert to units of £b)
    r_GDP = np.array(df.iloc[24,4:39])/100

    return date, sim_elms, r_GDP, num_years

def load_annual_GDP_ons():
    
    """
    Load in historical seasonally adjusted Gross Domestic Product (GDP)
    
    GDP data from the Office for National Statistics: https://www.ons.gov.uk/economy/grossdomesticproductgdp/timeseries/abmi/pn2

    Data is compiled using the gross fixed capital formation and is also used to estimate emissions from the other data sources.
    
    Inputs:
        
        None


    Outputs:
        
        - date, the year of the corresponding annual estimate. (type: datetime)
        
        - GDP, easonally adjusted Gross Domestic Product (GDP) (£billiion)
        
        - num_years, number of years of time-series. 
        

    """
    
    # Load in Data
    filename = 'series-170620.xls'
    excel_spreadsheet =  loaddata(filename,domain='economy')
    # sheet_names[0] only one sheet
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[0])
    # Get Year
    year = np.array(df.iloc[8:79,0])
    
    
    num_years = len(year)
    date = np.zeros(num_years,dtype=datetime.datetime)    

    for yr in range(0,num_years):
        
        date[yr] = datetime.datetime(int(year[yr]),1,1)  + (datetime.datetime(int(year[yr]),12,31) - datetime.datetime(int(year[yr]),1,1))/2.0
            
    # Time frame we are interested in
    elms = date > datetime.datetime(1990,1,1)

    # Get GDP (convert to units of £b)
    GDP = np.array(df.iloc[8:79,1])/1000

    return date[elms], GDP[elms], len(date[elms])

def load_precovid_emission_gov_estimate():
    
    """
    Load in future GHG emissions used by the UK goverment before COVID-19 outbreak. 
    
    Emission forecast from the Department for Business, Energy & Industrial Strategy: https://www.gov.uk/government/publications/updated-energy-and-emissions-projections-2018
    
    Inputs:
        
        None
    
    Outputs:
        
        - date, the year midpoint of the corresponding annual estimate. (type: datetime)
        
        - sim_elms, years that are simulated. (greater than 2016)
        
        - total_E_fore, the sum of GHG emissions across all the sectors. (units: MtCO2e)
        
        - total_E_95CI, 95% confidence interval of the forecast. (units: MtCO2e)
        
        - sector_E_fore, by sector GHG emissions. (units: MtCO2e)
        
        - sector_name, defined list of sectors.
        
        - num_years, number of years of time-series. 
    
    """
    
    # Load in Data
    filename = 'Copy of Annex-a-greenhouse-gas-emissions-by-source__16-May-2019_.xlsx'
    excel_spreadsheet =  loaddata(filename,domain='emissions')
    # sheet_names[1] provides for the reference scenarios. Breaks down GHG emission by sector (MtCO2e)
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[1])
    # Get Year
    year = np.array(df.iloc[10,1:])
    num_years = len(year)
    date = np.zeros(num_years,dtype=datetime.datetime)
    
    
    for yr in range(0,num_years):
                
        date[yr] = datetime.datetime(int(year[yr]),1,1)  + (datetime.datetime(int(year[yr]),12,31) - datetime.datetime(int(year[yr]),1,1))/2.0
            
    
    # Year of simulation    
    sim_elms = date>datetime.datetime(2017,1,1)
    
    # Defined Sectors
    sector_name = np.array(['Agriculture','Business','Energy supply','Industrial processes','LULUCF','Public','Residential','Transport','Waste management'])
    num_sector = len(sector_name)
    spreadsheet_rows = np.array(df.iloc[:,0])
    num_rows = len(spreadsheet_rows)       
    sector_E_fore = np.zeros((num_sector,num_years))
       
    
    # Find Emissions
    for sec in range(0,num_sector):
        
        for row in range(0,num_rows):
            
            if sector_name[sec] == spreadsheet_rows[row]:
                
                for yr in range(0,num_years):
                    
                    sector_E_fore[sec,yr] = np.array(df.iloc[row,1+yr])
                
                break
            
    # Find total emissions
    total_E_fore = np.sum(sector_E_fore,axis=0)
    # Get Uncertainty in forecast
    filename = 'Web_Figures_EEP2018.xlsx'
    excel_spreadsheet =  loaddata(filename,domain='emissions')
    # sheet_names[1] provides uncertainty
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[1])
    total_E_95CI = np.transpose(np.array(df.iloc[35:55,5:7]))
    a = total_E_fore[sim_elms]
    total_E_95CI[0,1] = a[0]
    total_E_95CI[1,1] = a[0]

    return date, sim_elms, total_E_fore, total_E_95CI, sector_E_fore, sector_name, num_years

def load_historical_emissions():
    
    
    """
    Load in the total historical GHG emissions for the UK.
    
    Emission data from the Department for Business, Energy & Industrial Strategy: https://www.gov.uk/government/statistics/final-uk-greenhouse-gas-emissions-national-statistics-1990-2017

    Inputs:
        
        None


    Outputs:
        
        - date, the year midpoint of the corresponding annual estimate. (type: datetime)
        
        - total_E_hist, the sum of GHG emissions across all the sectors. (units: MtCO2e)
        
        - sector_E_hist, by sector GHG emissions. (units: MtCO2e)
        
        - sector_name, defined list of sectors.
        
        - num_years, number of years of time-series. 
        
    
    """
    
    # Load in Data
    filename = 'Final_greenhouse_gas_emissions_tables_2017.xlsx'
    excel_spreadsheet =  loaddata(filename,domain='emissions')
    # sheet_names[3] breaks down GHG emission by sector (MtCO2e)
    sheet_names = excel_spreadsheet.sheet_names
    df = excel_spreadsheet.parse(sheet_names[3])
    # Get Year
    year = np.array(df.iloc[0,2:])
    num_years = len(year)
    date = np.zeros(num_years,dtype=datetime.datetime)
    for yr in range(0,num_years):
        
        date[yr] = datetime.datetime(int(year[yr]),1,1)  + (datetime.datetime(int(year[yr]),12,31) - datetime.datetime(int(year[yr]),1,1))/2.0
            
        
    # Defined Sectors
    sector_name = np.array(['Energy supply','Business','Transport','Public','Residential','Agriculture','Industrial processes','Land use, land use change and forestry','Waste management'])
    num_sector = len(sector_name)
    spreadsheet_rows = np.array(df.iloc[:,0])
    num_rows = len(spreadsheet_rows)       
    sector_E_hist = np.zeros((num_sector,num_years))

    # Find Emissions
    for sec in range(0,num_sector):
        
        for row in range(0,num_rows):
            
            if sector_name[sec] == spreadsheet_rows[row]:
                
                for yr in range(0,num_years):
                    
                    sector_E_hist[sec,yr] = np.array(df.iloc[row,2+yr])
                
                break
    # Find total emissions
    total_E_hist = np.sum(sector_E_hist,axis=0)
    
    return date, total_E_hist, sector_E_hist, sector_name, num_years


def savedata(df,filename='untitled',domain=None,directory=None,makedir=False):
    
    """
    Function saves panda dataframe in a csv format in a directory (default 
    assumes:‘…/code/data/output’). Domain corresponds to a given subfolder
    within directory.
    
    """
    
    
    if directory is None:
        
        directory = os.path.abspath(os.path.join(os.getcwd(), '..','data/output'))
  
    if domain is None:
        
        outputpath = directory 
        
    else:
          
        outputpath = os.path.abspath(os.path.join(directory,domain))

    if os.path.isdir(outputpath) == False:

        if makedir == False:

            print('Warning: "'+filename+'.csv" data has not been saved - directory "'+outputpath+'" does not exist.\n'+'Solution(s): either change argument makedir=True or find obtainable directory currently dataset to be saved in.')
        
            return
            
        else:
            
            os.makedirs(outputpath)
        
    if filename.find('csv') >= 0:
        
        outputfile = outputpath + '/' + filename
        
    else:
        
        outputfile = outputpath + '/'+filename+'.csv'
    
    df.to_csv(outputfile,index=False)
    
    return  


def loaddata(file,domain=None,directory=None):
    
    """
    Function returns either a panda dataframe or panda Excel type for 
    processing default assume this is based in local directory in 
    ‘…/code/data/raw’.
    
    """
    
    
    
    if directory is None:
    
        directory = os.path.abspath(os.path.join(os.getcwd(), '..','code/data/raw'))
        
    if domain is None:

        inputpath = directory 
        
    else:
          
        inputpath = os.path.abspath(os.path.join(directory,domain))

    if os.path.isdir(inputpath) == False:
        
        raise UserWarning('Error: directory "'+inputpath+'" does not exist.')

    inputfile = inputpath+'/'+file
    
    if os.path.exists(inputfile) == True:
        
        fmt = file.split('.')[-1]
        
    else:
        
        raise UserWarning('Error: file "'+file+'" does not exist in directory "'+inputpath+'".')
        
        
    panda_excelfmt = ['xls','xlsx','xlsm','xlsb','ods']
        
    if any(fmt in panda_fmt for panda_fmt in panda_excelfmt):
        
        if fmt == 'ods':
            
            return pd.read_excel(inputfile, engine="odf")
        
        else:
                
            return pd.ExcelFile(inputfile)
    
    else:
                
        raise UserWarning('Error: file format ".'+fmt+'" is unsupported by the function.\n'+'Solution(s): formats supported '+str(panda_excelfmt))
        


