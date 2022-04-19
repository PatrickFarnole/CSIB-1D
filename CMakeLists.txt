add_library(fabm_models_uvic OBJECT
            uvic_model_library.F90
            eco.F90
            icealgae.F90
            dic.F90
            dms.F90
            icedms.F90
           )

add_dependencies(fabm_models_uvic fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2016 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
#
           
            #npzd_papa.F90
            #dmsp_papa.F90
            