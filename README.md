# Analisis funcional de mutantes revertientes afectados en la ruta TOR de Magnaporthe oryzae
La cepa silvestre Guy11 de M.oryzae es sensible a la rapamicina y, en un medio que la contiene, exhibe un crecimiento significativamente reducido en comparación con un medio sin rapamicina. Sin embargo, el mutante Δrbp35 es capaz de crecer en presencia de rapamicina, demostrando resistencia al compuesto.  

Anteriormente a la realización de este trabajo, se habían seleccionado y secuenciado por DNA-seq un total de 6 estirpes revertientes procedentes del cultivo secuencial de Δrbp35 que habían mostrado diferencias morfológicas y colorimétricas y que revertieron la insensibilidad a la rapamicina. Se secuenciaron por DNA-seq también muestras de las estirpes control Guy11 y del parental mutante Δrbp35 tomadas en 3 tiempos diferentes para observar su evolución.

El mutante Δrbp35 exhibe alteraciones en la vía TOR y modificaciones en la expresión de transposones, pero la relación entre estos fenómenos no está completamente caracterizada. Dado que las cepas revertientes recuperaron la sensibilidad a la rapamicina, es razonable suponer que presentan alguna modificación genética vinculada a la vía TOR. Adicionalmente, se planteó la hipótesis de que estas cepas podrían tener un contenido de transposones diferente. 

Con el análisis de este experimento de DNA-seq se pretendían determinar las diferencias genómicas entre estirpes y relacionarlas con el cambio fenotípico observado. En concreto, se buscaba dilucidar si las regiones alteradas contenían genes relacionados con la ruta TOR y si los revertientes presentaban cambios en la cantidad de transposones. Además, este proyecto busca optimizar el flujo de trabajo seleccionando las herramientas bioinformáticas más adecuadas para el tratamiento de datos. 


### Orden de ejecución de los scripts
Durante este trabajo se ejecutaron los scripts:
1) analisis_calidad_inicial.sh
2) filtrado_trimmomatic.sh
3) filtrado_fastp.sh
4) mapeado_segemehl.sh
5) mapeado_bwa.sh
6) avg_ins_SOAPdenovo2_picard.sh
7) ensamblado_SOAPdenovo2.sh
8) ensamblado_SPAdes.sh
9) ensmablado_MaSuRCA.sh
10) ensamblado_MaSuRCA_Chromosome_scaffolder.sh
