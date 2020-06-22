# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 12:44:34 2020

@author: aargl
"""

import numpy as np
import datetime
from scipy.optimize import curve_fit


# Local Imports
from .modelling import *
from .secondary import *
from .ODE_model import run_simulations


class experiment:
    
    def __init__(self):
        
        """
        System class helps with the 
        
        """
        
        # Define temporal range
        self.timelimit = np.zeros(2,dtype=datetime.datetime)
        self.timelimit[0] = datetime.datetime(1990,1,1)
        self.timelimit[1] = datetime.datetime(2050,1,1)
        
        ### emissions        
        
        # Load in historical emissions by sector.
        self.year_E_hist_ccc, self.total_E_hist_ccc, self.sector_E_hist_ccc, self.sector_name_E, self.num_years_E_hist = load_historical_emissions()
        # Load in government emissions forecast.
        self.year_E_fore_gov, self.elms_E_fore_gov, self.total_E_fore_gov, self.total_E_fore_gov_95CI, sector_E_fore_gov, sector_name_E_fore_gov, self.num_years_E_fore = load_precovid_emission_gov_estimate()
                
        # Define Targets
        self.paris_wb_2 = (1-0.8) * self.total_E_hist_ccc[0]
        self.net_zero = 0 
        
        # Reorder to be consistent        
        self.num_E_sectors = len(self.sector_name_E)    
        org_hist_gov = np.zeros(self.num_E_sectors,dtype=int)        
        
        for sec1 in range(0,self.num_E_sectors):
            
            for sec2 in range(0,self.num_E_sectors):
                
                if self.sector_name_E[sec1] == sector_name_E_fore_gov[sec2]:
                    
                    org_hist_gov[sec1] = sec2
                    
                    break
                
                elif (self.sector_name_E[sec1] == 'Land use, land use change and forestry' and sector_name_E_fore_gov[sec2] == 'LULUCF'):
                    
                    org_hist_gov[sec1]  = sec2
                    
                    break
        
        
        # sector name short form (for plotting)
                
        self.sector_name_E_sf = np.array(['EneS','Bus','Tran','Pub','Res','Agr','IndP','LULUCF','WM'])
        self.sector_E_fore_gov = sector_E_fore_gov[org_hist_gov,:]

        ### economy
        
        # Load in historical economic data
        self.year_gdp_hist_ons, self.gdp_hist_ons, self.num_years_gdp_hist = load_annual_GDP_ons()
        self.month_gdp_hist_ons, self.gdp_hist_ons_mnth, self.num_months_gdp_hist = load_covid_GDP_ons() # By month
        
        # Load in projection of GDP
        self.year_gdp_fore_gov, self.elms_gdp_fore_gov, self.frac_change_gdp_fore_gove, self.num_years_gdp_fore = load_precovid_annual_GDP_estimate()
        
        # Find forecasted absolute values of GDP
        self.gdp_fore_gov = np.zeros(self.num_years_E_fore)
        
        for yr1 in range(0,self.num_years_E_fore):
            
            # Match Absolute Historical value to GDP before forecast year otherwise use percent change
            if self.elms_E_fore_gov[yr1] == False:
            
                for yr2 in range(0,self.num_years_gdp_hist):
                         
                    if self.year_gdp_hist_ons[yr1] == self.year_gdp_hist_ons[yr2]:
                        
                        self.gdp_fore_gov[yr1] = self.gdp_hist_ons[yr2]
                        
                        break
                    
            else:
                
                break
            
        self.gdp_fore_gov[self.elms_E_fore_gov] = self.gdp_fore_gov[self.elms_E_fore_gov==False][-1] * (1+np.cumsum(self.frac_change_gdp_fore_gove[self.elms_gdp_fore_gov]))
                    
        ###  carbon intensity (eta)
        self.eta_gdp_hist = np.zeros(self.num_years_E_hist)

        # Historical        
        for yr1 in range(0,self.num_years_gdp_hist):
            
            for yr2 in range(0,self.num_years_E_hist):
                
               if self.year_gdp_hist_ons[yr1] == self.year_E_hist_ccc[yr2]:
                    
                    self.eta_gdp_hist[yr1] = self.total_E_hist_ccc[yr2]/self.gdp_hist_ons[yr1]
                    
                    break
        # Forecast
        self.eta_gdp_fore = self.total_E_fore_gov/self.gdp_fore_gov
        
        ### COVID 19 Deaths
        self.day_mort_gov, self.mort_per_day_gov, self.R_per_day_gov = load_COVID19_mort()
        
        ### Google mobility
        self.day_m_google, self.m_google = load_google_mobility()
             
        self.tau_us, self.tau_as, self.Y_run, self.eta_run, self.E_run, self.R_run, self.m_run, self.mu_run = run_simulations()
        
        
        self.num_sec = len(self.tau_us)
        self.startdate_run = datetime.datetime(2020,1,31)
        self.endate_run = datetime.datetime(2050,1,1)
        self.num_days_run = len(self.Y_run[:,0])
        self.date_run = np.zeros(self.num_days_run,dtype=datetime.datetime)
        self.date_run[0] = self.startdate_run
        for day in range(1,self.num_days_run):
            
            self.date_run[day] = self.date_run[day-1] + datetime.timedelta(days=1)
        
        return
    
    
    def fit_lockdown_mobility(self):
        
        """
        Estimate the rate of change for mobility under lockdown
        
        Inputs:
            
        Outputs:
            
        
        
        """
        
        # Date of lockdown 
        self.date_lockdown = datetime.datetime(2020,3,26)
        # Assume public act before lockdown by two weeks
        self.date_public_lockdown = self.date_lockdown - datetime.timedelta(days=14)
        self.elms_m_lockdown = (np.gradient(self.m_week_50)<0) * (self.outbreak_obs_weekly>self.date_public_lockdown)
        self.outbreak_m_decrease_week = self.outbreak_obs_weekly[self.elms_m_lockdown]
        T = ((self.outbreak_m_decrease_week[-1] - self.outbreak_m_decrease_week[0]).days)/7+1

        self.t_lockdown_fit  = np.arange(0,T,step=1)/52
        
        popt, _ = curve_fit(model_lockdown_m,self.t_lockdown_fit,self.m_week_50[self.elms_m_lockdown],p0=[0.01,0.2],bounds=([0,0],[0.02,np.inf]))
    
        self.m_l = popt[0]
        self.tau_c = popt[1]
        
        
        return 
    
    
    def fit_easing_mobility(self):
        
        """
        Estimate the rate of change for mobility under lockdown
        
        Inputs:
            
        Outputs:
            
        
        
        """
        
        # Get start of increase in mobility

        self.outbreak_m_increase_week = self.outbreak_obs_weekly[self.elms_m_lockdown][-1]
        self.elms_m_easing = (self.outbreak_obs_weekly>=self.outbreak_m_increase_week)
        subset = self.outbreak_obs_weekly[self.elms_m_easing]
        T = ((subset[-1] - subset[0]).days)/7+1

        self.t_easing_fit  = np.arange(0,T,step=1)/52
        
        popt, _ = curve_fit(model_easing_m,self.t_easing_fit,self.m_week_50[self.elms_m_easing],p0=[1.0,1.0],bounds=([0,0],[np.inf,np.inf]))
     
        self.m_easing_0 = popt[0]
        self.tau_l = popt[1]  
        
        self.outbreak_m_increase_daily = self.outbreak_obs_days[self.outbreak_obs_days>=self.outbreak_m_increase_week]    
        T = ((subset[-1]-subset[0]).days)/365.25
        self.m_easing_fit = model_easing_m(T,self.m_easing_0,self.tau_l)
        print(self.m_easing_fit)
                
        
        return 
    
    def fit_covid_function(self):
        
        """
        Estimate the covid mobility rel
        """
        return
        
        
        
    
    def moving_avg_COVID19(self):
        
        """
        Estimate the week day average of COVID-19 data with confidence bounds.
        
        Inputs:
            
            None
            
        Outputs:
            
            - num_days_outbreak, the number of days with datapoints.
            
            - num_months_outbreak, the number of months with datapoints.
            
            - R_data_daily, R over the corresponding time series (shares the same temporal dimensions as m). 
            
            - m_data_daily, google mobility over the corresponding time series(shares the same temporal dimensions as R).
            
            - R_week_50, average of R the last week days.
            
            - R_week_95, 95% range of data used in week day avg R.
            
            - m_week_50, average of google mobility over the last 14 days.
            
            - m_week_95, 95% range of data used in 14 day avg m. (units: fraction change from baseline)
            
            - R_month_50, average of R in the month.
            
            - R_month_95, 95% confidence range of data used for month avg R.
            
            - m_month_50, average of google mobility in the month.
            
            - m_month_95, 95% confidence range of data used in month day avg m. (units: fraction change from baseline)
            
        """
        
        # Calculate moving weekly averages and range. Use the total range of the pandemic
        
        # First month of outbreak (that we have data for)
        first_month = min(self.day_mort_gov[0].month,self.day_m_google[0].month)
        # First day of outbreak
        first_day = min(self.day_mort_gov[0].day,self.day_m_google[0].day)
        self.outbreak_start_date = datetime.datetime(2020,first_month,first_day)
        # Last day of data
        last_data_month = max(self.day_mort_gov[-1].month,self.day_m_google[-1].month)
        last_data_day = max(self.day_mort_gov[-1].month,self.day_m_google[-1].month)
        self.outbreak_last_data_date = datetime.datetime(2020,last_data_day,last_data_day)
        
        self.num_days_outbreak = (self.outbreak_last_data_date-self.outbreak_start_date).days
     

        
        # Get days and data on days
        self.outbreak_obs_days = np.zeros(self.num_days_outbreak ,dtype=datetime.datetime)
        self.outbreak_obs_days[0] = self.outbreak_start_date
        self.t_outbreak = np.arange(0,self.num_days_outbreak,step=1)
        self.R_data_daily = np.nan*np.ones(self.num_days_outbreak )
        self.m_data_daily = np.nan*np.ones(self.num_days_outbreak )
        
        for day in range(0,self.num_days_outbreak):
            
            if day > 0:
                
                self.outbreak_obs_days[day] = self.outbreak_obs_days[day-1] + datetime.timedelta(days=1)
                
            for day2 in range(0,len(self.day_mort_gov)):
                
                if (self.outbreak_obs_days[day].day == self.day_mort_gov[day2].day and self.outbreak_obs_days[day].month == self.day_mort_gov[day2].month):

                    self.R_data_daily[day] = self.R_per_day_gov[day2]
                    
                
                    break
                
            for day3 in range(0,len(self.day_m_google)):
                
                if (self.outbreak_obs_days[day].day == self.day_m_google[day3].day and self.outbreak_obs_days[day].month == self.day_m_google[day3].month):
                    
                    self.m_data_daily[day] = self.m_google[day3]
                    
                    
                    break
                
                
        # Get weekly sets
                
        # Firstly we find weeks
        self.num_weeks_outbreak = 0

        for day in range(0,self.num_days_outbreak):

            if self.outbreak_obs_days[day].weekday() == 0:
                
                
                if day + 7 < self.num_days_outbreak-1:
                    
                    self.num_weeks_outbreak = self.num_weeks_outbreak + 1
                
        # Next find specific date for week
        self.outbreak_obs_weekly = np.zeros(self.num_weeks_outbreak,dtype=datetime.datetime)
        self.R_week_50 = np.nan*np.ones(self.num_weeks_outbreak)
        self.R_week_95 = np.nan*np.ones((2,self.num_weeks_outbreak))
        self.m_week_50 = np.nan*np.ones(self.num_weeks_outbreak)
        self.m_week_95 = np.nan*np.ones((2,self.num_weeks_outbreak))
        
        
        week = 0
        
        for day in range(0,self.num_days_outbreak):
                
            if self.outbreak_obs_days[day].weekday() == 0:
                
                
                if day + 7 < self.num_days_outbreak-1:
                    self.outbreak_obs_weekly[week] = self.outbreak_obs_days[day] + (self.outbreak_obs_days[day+7] - self.outbreak_obs_days[day])/2
                    self.R_week_95[0,week] = np.percentile(self.R_data_daily[day:day+8],5)
                    self.R_week_95[1,week] = np.percentile(self.R_data_daily[day:day+8],95) 
                    
                    self.R_week_50[week] = np.percentile(self.R_data_daily[day:day+8],50)
                    self.R_week_95[0,week] = np.percentile(self.R_data_daily[day:day+8],5)
                    self.R_week_95[1,week] = np.percentile(self.R_data_daily[day:day+8],95) 
                    self.m_week_95[0,week] = np.percentile(self.m_data_daily[day:day+8],5)
                    self.m_week_95[1,week] = np.percentile(self.m_data_daily[day:day+8],95)  
                    self.m_week_50[week] = np.percentile(self.m_data_daily[day:day+8],50)       
                    
                week = week + 1
            
            
                
        # Get Monthly sets
        # Firstly we find weeks
                
        self.num_months_outbreak = 0
        
        current_month = -1

        for day in range(0,self.num_days_outbreak):

            if self.outbreak_obs_days[day].month > current_month:
                
                current_month = self.outbreak_obs_days[day].month
                num_days_in_month = (datetime.datetime(2020,current_month+1,1))-datetime.datetime(2020,current_month,1)  
                self.num_months_outbreak = self.num_months_outbreak + 1
        
        
        # Next find specific date for week
        self.outbreak_obs_months = np.zeros(self.num_months_outbreak,dtype=datetime.datetime)
                                            
        self.R_month_50 = np.nan*np.ones(self.num_months_outbreak)
        self.R_month_95 = np.nan*np.ones((2,self.num_months_outbreak))
        self.m_month_50 = np.nan*np.ones(self.num_months_outbreak)
        self.m_month_95 = np.nan*np.ones((2,self.num_months_outbreak))
        
        
        current_month = -1
        month = 0
        
        for day in range(0,self.num_days_outbreak):
                
            if self.outbreak_obs_days[day].month > current_month:                
                    
                current_month = self.outbreak_obs_days[day].month
                dmonth = datetime.datetime(2020,current_month+1,1)-datetime.datetime(2020,current_month,1)
                self.outbreak_obs_months[month] = self.outbreak_obs_days[day] + (datetime.datetime(2020,current_month+1,1)-datetime.datetime(2020,current_month,1))/2
                num_days_in_month = min(day+dmonth.days,self.num_days_outbreak)
                self.R_month_95[0,month] = np.nanpercentile(self.R_data_daily[day:num_days_in_month],5)
                self.R_month_95[1,month] = np.nanpercentile(self.R_data_daily[day:num_days_in_month],95)
                self.R_month_50[month] = np.nanpercentile(self.R_data_daily[day:num_days_in_month],50)
                self.R_month_95[0,month] = np.nanpercentile(self.R_data_daily[day:num_days_in_month],5)
                self.R_month_95[1,month] = np.nanpercentile(self.R_data_daily[day:num_days_in_month],95) 
                self.m_month_95[0,month] = np.nanpercentile(self.m_data_daily[day:num_days_in_month],5)
                self.m_month_95[1,month] = np.nanpercentile(self.m_data_daily[day:num_days_in_month],95)  
                self.m_month_50[month] = np.nanpercentile(self.m_data_daily[day:num_days_in_month],50)       
                    
                month = month + 1
            
        return 
    
    def define_secdate(self):
        
        
        """
        Function defines the seculation date of model. We define a temporal resolution of 1 day from the midpoint of 2017 to 1/1/2050 and since 1990
        
        """
        
        # Since 2017
        self.start_date = datetime.datetime(2017,1,1) + (datetime.datetime(2017,12,31) - datetime.datetime(2017,1,1))/2 
        self.end_date = datetime.datetime(2050,1,1)
        self.ktime = (self.end_date - self.start_date).days + 1
        self.date = np.zeros(self.ktime,dtype=datetime.datetime)
        self.t = np.zeros(self.ktime)
        self.dt = 1/365.25
        
        for k in range(0,self.ktime):
            
            self.date[k] = self.start_date + datetime.timedelta(days=self.t[k]*365.25)

            if k < self.ktime-1:
                
                self.t[k+1] = self.t[k] + self.dt
                
        # Since 1990
        self.start_date_hist = datetime.datetime(1990,1,1) + (datetime.datetime(1990,12,31) - datetime.datetime(1990,1,1))/2     
        self.ktime_1990_2050 = (self.end_date - self.start_date_hist).days + 1
        self.date_1990_2050 = np.zeros(self.ktime_1990_2050,dtype=datetime.datetime)
        self.t_1990_2050 = np.zeros(self.ktime_1990_2050)
                
        for k in range(0,self.ktime_1990_2050):
            
            self.date_1990_2050[k] = self.start_date_hist + datetime.timedelta(days=self.t_1990_2050[k]*365.25)
            
            if (self.date_1990_2050[k].year == self.start_date.year and self.date_1990_2050[k].month == self.start_date.month and self.date_1990_2050[k].day == self.start_date.day):
                
                self.ktime_proj_crossing = k
        
                
            if k < self.ktime-1:
                
                self.t_1990_2050[k+1] = self.t_1990_2050[k] + self.dt    
        
        return
        
    def fit_gdp_noncovid_projection(self):
        
        """
        Function attempts to fit Gross Domestic Product (GDP) for the non-covid projection from 2017 and assume an exponential growth curve. 
        We currently assume zero volatility in this projection.
        
        Inputs:
            
            None
            
        Outputs:
            
            None
        
        """
        
        # Get First Year gdp of projection
        self.Y2017_noncovid_proj = self.gdp_fore_gov[self.elms_E_fore_gov][0]
        # Get temporal range in terms of years
        timedelta = self.year_E_fore_gov[self.elms_E_fore_gov] - self.year_E_fore_gov[self.elms_E_fore_gov][0]
        # Number of years over time
        num_years = len(timedelta)
        self.t_noncovid_gdpproj_fit = np.zeros(num_years)
        
        for yr in range(0,num_years):
            
            self.t_noncovid_gdpproj_fit[yr] = timedelta[yr].days/365.25
        
        
        popt, _ = curve_fit(model_expected_gdp,self.t_noncovid_gdpproj_fit,self.gdp_fore_gov[self.elms_E_fore_gov],p0=(self.Y2017_noncovid_proj,0.03),bounds=([self.Y2017_noncovid_proj*(0.999999),0.0],[self.Y2017_noncovid_proj*(1.0000001),0.1]))
        
        self.mu_noncovid = popt[1]
        
        # Project over fitted data
        self.Y_noncovid_proj_fit = model_expected_gdp(self.t_noncovid_gdpproj_fit,self.Y2017_noncovid_proj, self.mu_noncovid)
        
        # seculate over pro
        self.Y_noncovid = model_expected_gdp(self.t,self.Y2017_noncovid_proj,self.mu_noncovid)
        
        
        return
        
    def find_eta_projection(self):
        
        """
        Function fits the observed + projected trend in carbon intensity (eta)
        of GDP, outputs an assumed baseline intensity.
        
        Inputs:
            
            None
            
        Outputs:
            
            None
        
        """
        
        # Get temporal range in terms of years
        timedelta = self.year_E_fore_gov[self.elms_E_fore_gov] - self.year_E_fore_gov[self.elms_E_fore_gov][0]
        # Number of years over time
        num_years = len(timedelta)
        
        self.t_eta_fit = np.zeros(num_years)
        
        for yr in range(0,num_years):
            
            self.t_eta_fit[yr] = timedelta[yr].days/365.25
            
            
        popt, _ = curve_fit(model_expected_eta,self.t_eta_fit,self.eta_gdp_fore[self.elms_E_fore_gov],p0=(0.7,0.1,0.01))
        
        self.eta_0 = popt[0]
        self.eta_b = popt[1]
        self.xi = popt[2]
        self.eta = model_expected_eta(self.t,self.eta_0,self.eta_b,self.xi)
        
        self.E_noncovid = model_emissions(self.eta,self.Y_noncovid)
        
        return
        

        
                
        
        
        
        
     





#%%
    #%% Plotting 
        

        
    
