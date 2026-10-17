# PyWKSpectra

<p align="center">
  <img src="assets/pywkspectra_logo.png" alt="PyWKSpectra Logo" width="280">
</p>

**Wheeler–Kiladis Space–Time Spectral Analysis in Python**

*A focused Python package for investigating tropical atmospheric variability through wavenumber–frequency spectra.*

## About PyWKSpectra
Aquí explicar un poco lo que hace este paquete y su filosofía

## Installation and examples

Explicar como instlarlo usando pip (supon que ya está disponible) y que paquetes se deben instalar. 
Agrega tambien las que permitan usar grib 

Ademas menciona el ejemplo y que usa los datos de ejemplos e indicr que paquetes se deben tener tambien para poder correr el cuadernos de jupyter.



## A little history (and a note on attribution)

This project began long before the name *PyWKSpectra* existed.

Este trabajo comenzó durante mi doctorado en 2014 cuando me encontraba en una estancia de investigación con el maravilloso George Kiladis del PSD de la NOAA. La idea original de este trabajo era construir una implemenciación en python del código disponible en NCL. The original work was inspired by the **NCAR Command Language (NCL)** implementation of Wheeler–Kiladis spectral analysis, especially its [`wkSpaceTime` diagnostic](https://www.ncl.ucar.edu/Document/Functions/Diagnostics/wkSpaceTime.shtml). That implementation was an important reference for bringing the approach to Python.

Sin embargo, el proyecto evolucionó y paso de depender menos de la forma en que NCL hacía el cálculo a usar formas mas pythonisticas. El codigo evolucionó en un paquete llamado wk_spectra, sin embargo la altisima carga acadmica y de trabajo y el nulo appoyo insittucional no me perm itió avanzar en desarrollar una version mejor documentada. El paquete estaba llamando la atención pero lamentablemente no tenia tiempoo disponible para poder estar atendiendo las oclicitudes y explicando el funcionameinto del paquete. Ademas me desilucionó ver que habían muchos usos del paquete pero poco reconocimiento a las cientos de horas que investí en el desarrolo. Lamentablemewnte el parque wk_spectra ha sido usado por muchos trabajos pero en las citaciones han atribuido el paquete a otro autor.

Aqui podemos agregar cosas de 

I subsequently developed **[`wk_spectra`](https://github.com/mmaiergerber/wk_spectra)**, an earlier Python implementation of the method. Its original authorship and copyright records identify **Alejandro Jaramillo** as the original developer. I appreciate the interest and contributions that helped the software reach other researchers.

Over time, maintaining and extending the package without sustained support became difficult. Its development history also became less visible than I had hoped, sometimes leading to confusion about the original authorship. This is not meant as criticism of any individual; it is simply a reminder that **software provenance deserves the same care as data provenance**.

After all, Git can track every changed line of code, but citations still need a little help from humans.

**PyWKSpectra** is a new effort to improve the scientific implementation, modernize the software, simplify its use, and make the history, authorship, and citation instructions explicit. Previous contributions will be acknowledged appropriately.

Esto que sigue tambien puede ser integrado. 
Tambien decir que me pueden contactar con preguntas y sugerencias pero les pido disculpas si tardo en responder o si no respkindo. Atenderé solo solicitudes rasonables ya que mis cargas academicas, y adminisrativas son demasiado elevadas y todo esto lo hago en mi poco tiempo libre. 



## How to cite PyWKSpectra

**If PyWKSpectra contributes to your figures, analyses, or scientific results, please cite the software and its author.**

The foundational Wheeler and Kiladis (1999) article explains the scientific method; a software citation additionally recognizes the implementation, verification, documentation, and maintenance needed to make that method accessible.

A formal `CITATION.cff` file and a DOI for an archived software release will be provided before the first public release. Until then, please acknowledge the package as **PyWKSpectra — Alejandro Jaramillo Moreno** and link to its repository where appropriate.

*Spectral power may be symmetric or antisymmetric. Authorship should not be ambiguous.*

## Publications using PyWKSpectra

**Have you used PyWKSpectra in a paper, preprint, thesis, or report? Please send me the reference!**

I would be delighted to read your work and include it in a growing list of publications using the package. Seeing what researchers discover with the software is one of the most rewarding parts of sharing it.

*A publication list will be added once references are available.*

## The best way to say thanks

PyWKSpectra is being developed to support scientific research. If it saves you time or contributes to a result, there are three especially welcome ways to say thanks:

1. **Cite PyWKSpectra** in your publication.
2. **Share your article** so I can add it to the list of research using the package.
3. **Invite me to collaborate** when there is a meaningful scientific opportunity to work together.

Coffee certainly helps keep scientific code running, but **citations and collaborations are better long-term fuel**.

Constructive feedback, bug reports, and appropriately credited contributions are welcome too.

## Scientific reference

Wheeler, M., and G. N. Kiladis (1999): Convectively coupled equatorial waves: Analysis of clouds and temperature in the wavenumber–frequency domain. *Journal of the Atmospheric Sciences*, **56**, 374–399.

NCL reference: [`wkSpaceTime`](https://www.ncl.ucar.edu/Document/Functions/Diagnostics/wkSpaceTime.shtml).

## Author and maintainer

**Alejandro Jaramillo Moreno**  
Departamento de Cambio Climático  
Instituto de Ciencias de la Atmósfera y Cambio Climático (ICAyCC)  
Universidad Nacional Autónoma de México (UNAM)  
Mexico City, Mexico

**ORCID:** [0000-0002-4175-7726](https://orcid.org/0000-0002-4175-7726)  
**Email:** [ajaramillo@atmosfera.unam.mx](mailto:ajaramillo@atmosfera.unam.mx)

## License and acknowledgments

A distribution license will be selected before the first public release. Until that is done, **no general permission to redistribute or reuse new PyWKSpectra code is granted by this README**. Existing code reused from earlier projects remains subject to its original license terms and attribution requirements.

PyWKSpectra acknowledges the pioneering work of **Wheeler and Kiladis (1999)**, the **NCL** implementation that inspired the original Python work, and the contributions of researchers and users over the years.

---

*May your spectral peaks be physically meaningful, your tests pass, and your citations point to the right author.*
