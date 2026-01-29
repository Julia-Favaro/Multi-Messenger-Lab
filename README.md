# Julia Favaro's Project - Multi Messenger Laboratory 2023-2024

This course was designed by Prof. Massimiliano Razzano and Prof. Barbara Patricelli for the Master's degree in Physics at the University of Pisa.
University of Pisa webpage: [https://www.df.unipi.it/](https://www.df.unipi.it/)

Syllabus of this class: [Astrophysics and Multimessenger Laboratory](https://unipi.coursecatalogue.cineca.it/corsi/2024/10452/insegnamenti/2025/52570_695986_78712/2023/52570?annoOrdinamento=2023)

### Repository structure
- _Prel_ folder: preliminary lab experience on OOP, Pandas library and FITS files (individual work) .
- _Radio_ folder: Lab experience on radio mapping the Milky Way with radiotelescopes operating at 21 cm and creating the rotational velocity curve of the galaxy.
- _Optical_ folder: Lab experience on Photometry and CCD data reduction using the Lab Telescopes and ATIK 460 CCD
- _GrB_ folder:  Lab experience verifying the light curve of a GRB event using the FERMI GRB dataset. 
- _GWs_ folder: Lab experience on data visualization of a GW event through LIGO Hanford LIGO Livingston and Virgo datasets.(individual work)

Each folder contains a README file that goes through all the code files and the results obtained.

### Objectives of the course
- Understand the main experimental techniques for detecting electromagnetic signals and gravitational waves;
- Understand how to access and use public astrophysical data archives;
- Perform observations and data reduction operations with real instruments and datasets;
- Analyze and interpret data from radio, optical, X-ray, gamma-ray, and gravitational-wave observations;
- Design and implement data analysis activities using tools and methodologies employed in contemporary multimessenger astrophysical research.

### Usage, data, and libraries
My partners and I obtained the datasets during our lab lectures. Documentation is embedded in the code. 

**Prerequisites**: 

-   Python 3.9 or higher.
-   A virtual environment manager like Conda is recommended.
-   For the GRB experience, a different approach was used using the [HEASARC](https://heasarc.gsfc.nasa.gov/docs/software/index.html) tools provided by NASA Fermi-BGM 

**Dependencies**:

    - numpy - array and math functions
    - scipy - numerical methods 
    - matplotlib- plotting and data visualization
    - pandas - for working with datasets structures
    - astropy - astronomy utilities (physical units, coordinates, time-series, FITS I/O)
    - photoutils -  photometric analysis of astronomical images
    - GWpy -analysis and visualization of gravitational‑wave data (LIGO/Virgo/KAGRA)

The analysis was executed in Jupyter notebooks located in each directory.

### Authors and acknowledgment
- Julia Favaro's bio: I'm a master's student from the University of Pisa (Italy). I am currently pursuing a career in high energy astrophysics and data analysis.
- Lab partners: Francesco Certaldi and Chiara Masia, master's students from the University of Pisa.
- Supervisors and professors: Prof. Massimiliano Razzano and Prof. Barbara Patricelli.

### Support
For more explanation or any suggestions, write me at j.favaro@studenti.unipi.it .

### License
GNU General Public License v3.0 or later

### Project status
Final exam in July 2024. Final grade: 30/30. 
