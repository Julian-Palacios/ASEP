# ASEP
Repositorio para realizar Análisis Sísmico de Edificaciones en Python

## ¿Qué realiza los códigos?
Realiza análisis estático y dináminco modal espectral de una edificación de 'n' niveles según la normativa sismorresistente (NTP E030). Estos análisis se realizan a través de la libreria OpenSeesPy y otras funciones programadas en este repositorio.
Además, crea imágenes y genera un reporte automatizado en word.

## Requisitos:
- Instalar openssespy (pip install openseespy)
- Para ejecutar adecuadamente el código, es necesario comentar en la función plot_modeshape de la librería openseespy, la línea que plotea *plt.show()*.
- Instalar python-docx (pip install python-docx)
- Adicionalmente, es necesario usar librerías como *matplolib*, *numpy* y *pandas*.

## Consultas:
- jpalaciose@uni.pe
