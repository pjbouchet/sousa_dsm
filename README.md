## Contrasting aerial abundance estimates for two sympatric dolphin species at a regional scale using distance sampling and density surface modelling

This repository contains R code accompanying the article:

Raudino HC, Bouchet PJ, Douglas C, Douglas R, Waples K. Contrasting aerial abundance estimates for two sympatric dolphin species at a regional scale using distance sampling and density surface modelling.

#### Abstract

Monitoring wildlife populations over scales relevant to management is critical to supporting conservation decision-making in the face of data deficiency, particularly for rare species occurring across large geographic ranges. The Pilbara region of Western Australia is home to two sympatric and morphologically similar species of coastal dolphins --- the Indo-pacific bottlenose dolphin (Tursiops aduncus) and Australian humpback dolphin (Sousa sahulensis) --- both of which are little known, believed to be declining in numbers, and facing increasing pressures from the combined impacts of environmental change and rampant industry activities. The aim of this study was to develop spatially explicit models of bottlenose and humpback dolphin abundance in Pilbara waters that could inform decisions about coastal development at a regional scale. Aerial line transect surveys were flown from a fixed-wing aircraft in the austral winters of 2015, 2016, and 2017 across a total area of 33,420 km<sup>2</sup>. Spatio-temporal patterns in dolphin density were quantified using a density surface modelling (DSM) approach, accounting for imperfect detection as well as both perception and availability bias. We estimated the abundance of bottlenose dolphins at 3,713 (95% CI = 2,679---5,146; average density of 0.189 ± 0.046 SD individuals per km<sup>2</sup>) in 2015, 2,638 (95% CI = 1,670---4,168; 0.159 ± 0.135 individuals per km<sup>2</sup>) in 2016 and 1,635 (95% CI = 1,031---2,593; 0.101 ± 0.103 individuals per km<sup>2</sup>) in 2017. No humpback dolphins were detected in 2015, but their estimated abundance was 1,546 (95% CI = 942---2,537; 0.097 ± 0.03 individuals per km<sup>2</sup>) and 2,690 (95% CI = 1,792---4,038; 0.169 ± 0.064 individuals per km<sup>2</sup>) in 2016 and 2017, respectively. Dolphin densities were greatest along the inshore belt, with hotspots in Exmouth Gulf, the Dampier Archipelago, and Great Sandy Islands. Our results provide a benchmark on which future risk assessments can be based to better understand the overlap between anthropogenic pressures and important dolphin habitats in tropical northwestern Australia.

#### Funding & Acknowledgements

We are grateful to the pilots that kept us safely on transect and got us back down on the ground included Simon Garland, Steve Ewen and Gail Nolan flying for Shine Aviation Services. Kym Reeve successfully coordinated all three aerial surveys and ensured the high quality of the data collected in the field. Sincere thanks to the observers Christophe Cleguer, Daniella Hanf, Jane Kennedy, Krista Nicholson, Margie Rule, Chandra Salgado-Kent, Jenny Smith, Julian Tyne, Kristel Wenziker and Erin Wyatt. This research was funded by the Chevron-operated Wheatstone LNG Project's State Environmental Offsets Program and administered by the Department of Biodiversity, Conservation and Attractions. The Wheatstone Project is a joint venture between Australian subsidiaries of Chevron, Kuwaut Foreign Petroleum Exploration Company (KUFPEC), Apache Corporation and Kyushu Electric Power Company, together with PE Wheatstone Pty Ltd (part-owned by TEPCO). This project was also supported by Woodside through the Pluto LNG Environmental Offsets Program.

#### Files

`sousa_dsm_main.R`: Main script for the DSM analysis

`sousa_dsm_functions.R`: Supporting functions

`sousa_dsm_figures.R`: Code used to generate figures for the paper

`sousa_dsm_mrds.R`: Perception bias estimation

`sousa_dsm_truncation.R`: Code used to estimate right-truncation distances matching the time windows used in availability factors.
