# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 12:33:41 2020

@author: aargl
"""

from subroutines.primary import *
from subroutines.plotting import *


if __name__ == '__main__':
    
    """
    Runs the COVID / Economy / Emissions Experiment
    
    """

    run = experiment()
    run.define_secdate()
    run.fit_gdp_noncovid_projection()
    run.find_eta_projection()
    run.moving_avg_COVID19()
    run.fit_lockdown_mobility()
    run.fit_easing_mobility()
    run.fit_covid_function()
    #plot_fig1(run)
    plot_fig2(run)
    plot_fig3(run)

