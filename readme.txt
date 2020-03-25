{\rtf1\ansi\ansicpg1252\cocoartf2511
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red34\green139\blue34;}
{\*\expandedcolortbl;;\csgenericrgb\c13333\c54510\c13333;}
\paperw11900\paperh16840\margl1440\margr1440\vieww14660\viewh8840\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs20 \cf2 \

\fs24 Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems\cf0 \
\
\cf2 Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,\cf0  \cf2 Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.\cf0 \
\cf2 (2020) Proceedings of the National Academy of Sciences\cf0 \
\cf2  \cf0 \
\cf2 Barcelona,Spain. March, 2020.\cf0 \
\cf2  \cf0 \
\
\
First, run LEiDA_PsiloData.m to extract the probabilistic metastable substate (PMS) space for the empirical psilocybin BOLD data (placebo and active condition) using the leading eigenvector dynamics analysis (LEiDA) method.\
\
Second, run optim_placebo_psilo.m for fitting the whole-brain model to the PMS space of the placebo condition of psilocybin using only the neuronal system. Use the output from LEiDA_PsiloData.m. This yields a global neuronal coupling parameter (G=1.6). Use read_G_cluster.mat to visualize the results.\
\
Third, proceed to study the role of the neurotransmitter dynamical system by coupling both systems in the model run neuromodulador2D_psilo_raphenorm.m using the output from LEiDA_PsiloData.m and incorporating the connectivity of the raphe nucleus. Use read_Gw_2D_cluster.m to visualize the results.\
\
Further, for finding the optimum fit of mutually coupled whole-brain model as a function of excitatory and inhibitory coupling parameters for when the system is disconnected (Ws E and I =0), use psilo_w_0.m, and use psilo_w_optimum.m for finding the most significantly optimal fit to the empirical PMS. Use read_single_w.m to visualize both outputs.\
\
Finally, for the various manipulations to the model use manipulation_psilo_raphenorm.m and read_manipulation.m for visualization.\
\
 }