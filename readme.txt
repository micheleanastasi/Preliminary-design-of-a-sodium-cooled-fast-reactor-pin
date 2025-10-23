++ LOOK AT THE PDF FILE IF INTERESTED ONLY IN FINAL RESULTS (JUST A 5-PAGE SUMMARY AS REQUIRED), OTHERWISE LOOK AT CODE+PLOTS! ++

PRELIMINARY DESIGN OF A SODIUM-COOLED FAST REACTOR PIN
That was a project group which I worked for, contribuiting heavily on the design/coding part.

GUIDE:
Main script:
- main_th.py --> main script for thermal, remember to switch LoadExisting according to needs
- main_mech.py --> main script for mechanical, remember to switch LoadExisting according to needs
- Thermo-mechanics_preliminary.ipynb --> preliminary sizing (see report); NOTE: open it with COLAB


- excels --> some excels from testmain_hot_geometry.py depicting various results and properties
- figures --> important! Results from main scripts depicting results
- functions --> folder where there are important functions used to compute everything
- main_xxx_saves --> folders containing saved results as .npy, to be opened only by the code
- restr_saves --> columnar and void radius contained here, generated at low burn-up and used for higher ones

- testmain_hot_thickness_margin_calc.py --> used to compute optimum for thickness vs margin to fuel melting
- testmain_xxx.py --> various script (not the final result) used to perform calcs and tests
