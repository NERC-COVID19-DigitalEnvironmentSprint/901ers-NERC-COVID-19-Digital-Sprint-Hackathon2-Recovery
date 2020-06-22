# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 16:18:47 2020

@author: -
"""



import os
from matplotlib.ticker import MultipleLocator
import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from matplotlib import pyplot as plt
from matplotlib import gridspec
import datetime

# Local imports
from .primary import experiment

    
def plot_fig1(model,filename=None,output_dir=None,fmt=None):
    
    if output_dir is None:
        
        output_dir =  os.path.abspath(os.path.join(os.getcwd(), '..','presentation/'))
        
    if filename is None:
        
        filename = 'fig1'
        
    if fmt is None:
        
        fmt = 'png'
        
        
    outputfile = output_dir + filename +'.'+ fmt

    Month = mdates.MonthLocator()   
    xfmt = mdates.DateFormatter('%b')
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(3,1,height_ratios=[1,1,1],hspace=0.2)
    ax = fig.add_subplot(gs[0,0])    
    ax.scatter(model.outbreak_obs_days,model.m_data_daily,marker='.',color=[0,0,0,0.25],label='Google Daily')
    ax.scatter(model.outbreak_obs_weekly[model.elms_m_easing],model.m_week_50[model.elms_m_easing],marker='^',edgecolor=[0,0,0,1],facecolor=[0,0,0,0],s=60,label='Weekly Median')
    ax.scatter(model.outbreak_obs_weekly[model.elms_m_lockdown],model.m_week_50[model.elms_m_lockdown],marker='v',edgecolor=[0,0,0,1],facecolor=[0,0,0,0],s=60)
    ax.scatter(model.outbreak_obs_weekly[~(model.elms_m_lockdown+model.elms_m_easing)],model.m_week_50[~(model.elms_m_lockdown+model.elms_m_easing)],marker='o',edgecolor=[0,0,0,1],facecolor=[0,0,0,0],s=60)    
    ax.plot(model.outbreak_m_increase_daily,model.m_easing_fit)
    ax.plot([model.date_public_lockdown,model.date_public_lockdown],[0,1.2],color='#7570b3',linewidth=1.5)
    ax.annotate("`Public' Lockdown",xy=(model.date_public_lockdown-datetime.timedelta(days=23),0.4),color='#7570b3')
    ax.plot([model.date_lockdown,model.date_lockdown],[0,1.2],color='#d95f02',linewidth=1.5)
    ax.annotate('UK Gov Lockdown',xy=(model.date_lockdown+datetime.timedelta(days=1),1),color='#d95f02')
    ax.plot([model.date_lockdown,model.date_lockdown],[0,1.2],color='#d95f02',linewidth=1.5)
    ax.plot([model.outbreak_m_increase_week,model.outbreak_m_increase_week],[0,1.2],color='#1b9e77',linewidth=1.5)
    ax.annotate('Increase Mobility',xy=(model.outbreak_m_increase_week+datetime.timedelta(days=1),0.8),color='#1b9e77')
    ax.xaxis.set_major_locator(Month)
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xticks(model.outbreak_obs_weekly,minor=True)
    #ax.yaxis.set_minor_locator(MultipleLocator(50))
    #ax.set_xlim(focused_startdate,focused_endate)
    ax.set_ylim(0.2,1.15)
    ax.set_ylabel(r'$m$ (Baseline $m=1$)')
    ax.set_title(r'(a) UK Mobility')    
    ax.grid()
    
    plt.savefig(outputfile)
#        self.outbreak_m_increase_week = self.outbreak_obs_weekly[self.elms_m_lockdown][-1]
#        self.elms_m_easing = (self.outbreak_obs_weekly>=self.outbreak_m_increase_week)
#        subset = self.outbreak_obs_weekly[self.elms_m_easing]
#        T = ((subset[-1] - subset[0]).days)/7+1
#
#        self.t_easing_fit  = np.arange(0,T,step=1)/52
#        
#        popt, _ = curve_fit(model_easing_m,self.t_easing_fit,self.m_week_50[self.elms_m_easing],p0=[1.0,1.0],bounds=([0,0],[np.inf,np.inf]))
#    
#        
#        self.m_easing_0 = popt[0]
#        self.tau_l = popt[1]  
    return
        
    
    
def plot_fig2(model,output_dir=None,filename=None,fmt=None):
    
    if output_dir is None:
        
        output_dir =  os.path.abspath(os.path.join(os.getcwd(), '..','code/presentation/'))
        
    if filename is None:
        
        filename = 'fig2'
        
    if fmt is None:
        
        fmt = 'png'
        
        
    outputfile = output_dir + filename +'.'+ fmt
    
    
        
    years = mdates.YearLocator()   
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(2,2,height_ratios=[1,1],width_ratios=[1,1],hspace=0.2,wspace=0.2)
    ax = fig.add_subplot(gs[0,0])    
    
    for sec in range(0,model.num_sec):
        
        ax.scatter(model.tau_us[sec],model.tau_as[sec])

    ax.set_xlim(0,3.5)
    ax.set_ylim(0,3.5)
    ax.set_ylabel('Rate of Adaptation of Economy (%/yr)')
    ax.set_xlabel('Rate of Easing of Lock Down (%/yr)')
    ax.set_title('(a) Picked Scenarios')
    ax.grid(axis='both',which='major',zorder=0)
    
    ax = fig.add_subplot(gs[0,1])
    focused_endate = datetime.datetime(2025,1,1) 
    focused_startdate = datetime.datetime(2020,1,1)
    for sec in range(0,model.num_sec):
        
        ax.plot(model.date_run,model.R_run)
        
    ax.set_ylim(0,4)
    ax.set_ylabel('R')
    ax.set_title('(b) Rate of Reproduction')
    ax.xaxis.set_minor_locator(years)
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.grid(axis='both',which='major',zorder=0)
    ax.set_xlim(focused_startdate,focused_endate)
    
    ax = fig.add_subplot(gs[1,0])
       
    for sec in range(0,model.num_sec):
        
        ax.plot(model.date_run,model.mu_run*100)
    
    ax.xaxis.set_minor_locator(years)
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.set_xlim(focused_startdate,focused_endate)
    ax.set_ylabel(r'$\mu$ (%/yr)')
    ax.set_title(r'(c) GDP Growth Rate')
    ax.grid()
    ax = fig.add_subplot(gs[1,1])
    
    for sec in range(0,model.num_sec):
        
        ax.plot(model.date_run,model.m_run)
    
    ax.xaxis.set_minor_locator(years)
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.set_xlim(focused_startdate,focused_endate)
    ax.set_ylim(0,1)
    ax.set_ylabel(r'$m$ (Baseline $m=1$)')
    ax.set_title(r'(d) Mobility')    
    ax.grid()
    plt.savefig(outputfile)
    return

def plot_fig3(model,output_dir=None,filename=None,fmt=None):
    
    if output_dir is None:
        
        output_dir = os.path.abspath(os.path.join(os.getcwd(), '..','code/presentation/'))
        
    if filename is None:
        
        filename = 'fig3'
        
    if fmt is None:
        
        fmt = 'png'
        
        
    outputfile = output_dir + filename +'.'+ fmt
    
    years = mdates.YearLocator()   
    fig = plt.figure(figsize=(12,12))
    gs = gridspec.GridSpec(3,2,height_ratios=[1,1,1],width_ratios=[1,0.2],hspace=0.2)
    
    ax = fig.add_subplot(gs[0,0])    
    ax.plot(model.year_gdp_hist_ons,model.gdp_hist_ons,color='black',label='ONS: Seasonally adjusted GDP (CV measure)',zorder=3,linewidth=2)
    ax.plot(model.year_E_fore_gov[model.elms_E_fore_gov],model.gdp_fore_gov[model.elms_E_fore_gov],label='UK Gov Proj 2018',color='#808080',zorder=3,linewidth=2)
    ax.plot(model.month_gdp_hist_ons,model.gdp_hist_ons_mnth,color='red',linewidth=1.5,label='ONS: Monthly Measure',zorder=3)
    ax.plot(model.date,model.Y_noncovid,linewidth=1.5,color='blue',label='Modelled Pre-Covid 19 Growth ('+str(round(model.mu_noncovid*100,2))+r'$\ \mathrm{\%\ yr^{-1}}$)',zorder=2)
    
    for sec in range(0,model.num_sec):
        
        ax.plot(model.date_run,model.Y_run[:,sec])
    
    ax.set_xlim(model.timelimit)
    ax.set_ylabel('£ billion')
    ax.set_title('(a) UK GDP')
    ax.xaxis.set_minor_locator(years)
    ax.set_ylim(0,4000)
    ax.set_yticks([0,500,1000,1500,2000,2500,3000,3500,4000])
    ax.yaxis.set_minor_locator(MultipleLocator(250))
    ax.legend(fontsize=8,framealpha=1)
    ax.grid(axis='both',which='major',zorder=0)
    ax.xaxis.set_major_formatter(NullFormatter())
    
    
    
    ax = fig.add_subplot(gs[1,0])    
    ax.plot(model.year_E_hist_ccc,model.total_E_hist_ccc,label='(c) UK GHG Nat.Statistics: 1990-2017',color=[0,0,0],zorder=3)
    elms = model.elms_E_fore_gov
    elms[26] = True
    ax.plot(model.year_E_fore_gov[elms],model.total_E_fore_gov[elms],label='UK Gov Proj 2018 (Ref Case with 95% CI)',color='#808080',zorder=3)
    ax.plot(model.year_E_fore_gov[model.elms_E_fore_gov],model.total_E_fore_gov_95CI[0,:],ls='-.',color='#808080',zorder=2)
    ax.plot(model.year_E_fore_gov[model.elms_E_fore_gov],model.total_E_fore_gov_95CI[1,:],ls='-.',color='#808080',zorder=2)   
    ax.plot(model.date,model.E_noncovid,linewidth=1.5,color='blue',label='Modelled Pre-Covid 19 Emissions')
    ax.plot(model.timelimit,[model.paris_wb_2,model.paris_wb_2],color='#7570b3',zorder=2)
    
    for sec in range(0,model.num_sec):
        
        ax.plot(model.date_run,model.E_run[:,sec])
        
    ax.annotate(r'Paris $2^\degree$C Target',xy=(datetime.datetime(2000,1,1),model.paris_wb_2*1.1),fontsize=8,color='#7570b3',zorder=2)
    ax.annotate('Net Zero',xy=(datetime.datetime(2000,1,1),model.paris_wb_2*1.1-model.paris_wb_2),fontsize=8,color='#1b9e77',zorder=2)
    ax.plot(model.timelimit,[model.net_zero,model.net_zero],color='#1b9e77',zorder=2)
    ax.set_xlim(model.timelimit)
    ax.set_ylim(-5,810)
    ax.set_ylabel('MtCO2e')
    ax.set_title('(c) UK GHG Emissions')
    ax.legend(fontsize=8,framealpha=1)
    ax.xaxis.set_minor_locator(years)
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.grid(axis='both',which='major',zorder=0)
    
    
    
    ax = fig.add_subplot(gs[2,0])
    elms = model.elms_E_fore_gov
    elms[26] = True
    ax.plot(model.year_E_hist_ccc,model.eta_gdp_hist,color='black',label='Historical Emissions:GDP ratio')
    ax.plot(model.year_E_fore_gov[model.elms_E_fore_gov],model.eta_gdp_fore[model.elms_E_fore_gov],color='#808080',label='Projected Emissions:GDP ratio',zorder=3)
    ax.plot(model.date,model.eta,linewidth=1.5,color='blue',label='Modelled Carbon Intensity')

    ax.set_xlim(model.timelimit)
    ax.set_ylabel('MtCO2e / £ billion')
    ax.set_title('(d) UK Carbon Intensity')
    ax.legend(fontsize=8,framealpha=1)
    ax.xaxis.set_minor_locator(years)
    ax.grid(axis='both',which='major',zorder=0)
    ax.set_ylim([-0.1,0.71])  
    
    plt.savefig(outputfile,dpi=500)
    
    

    return