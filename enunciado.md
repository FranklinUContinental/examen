# Examen de Primera Unidad — Enunciado

**Curso:** Métodos Estadísticos para Investigaciones Experimentales
**Programa:** Escuela de Posgrado — Unidad de Posgrado en Ciencias Físicas y Matemáticas
**Universidad:** Universidad Nacional de Trujillo
**Docente:** Prof. Roger Reyna Segura

---

## 1. Consigna del examen

Se entrega al estudiante la tesis titulada:

> **«Radiación UV-C para inactivar patógenos alimentarios *Escherichia coli* y *Salmonella typhimurium* en carne de cuy contaminado»**

cuyo problema central es determinar **de qué manera la radiación UV-C influye en la inactivación de los patógenos alimentarios *E. coli* y *S. typhimurium*** inoculados en carne de cuy (*Cavia porcellus*).

A partir de esa tesis, el examen pide (textualmente, de la presentación del docente):

1. **Evaluar si el diseño y los resultados de la investigación son correctos.**
2. **Determinar qué diseño estudiado en la primera unidad sería el mejor** para esta investigación, **planteando el modelo y las hipótesis** correspondientes.
3. **Evaluar los supuestos básicos** con el modelo seleccionado (considerando cada paso).
4. **Realizar un nuevo análisis estadístico**, describiendo adecuadamente los resultados.
5. **Formular las conclusiones finales** de la investigación.

> En síntesis, no se trata de "repetir" la tesis, sino de **auditarla estadísticamente**: criticar el diseño declarado, proponer el modelo correcto, contrastar hipótesis, validar supuestos, rehacer el ANOVA y concluir.

---

## 2. Contexto de la investigación original (tesis)

La tesis aplica un sistema de irradiación UV-C como **tecnología emergente de conservación** de carne de cuy, alimento de alto valor proteico (≈19.5 % de proteína) y creciente interés comercial en Perú.

**Montaje experimental descrito en la tesis:**

- Cámara de irradiación de **50 × 30 × 50 cm**, recubierta internamente con papel de aluminio para maximizar el aprovechamiento de la radiación.
- Lámpara **UV-C** montada en la parte superior, con intensidad de radiación **0.654 ± 0.04 J/cm²**.
- Carne de cuy cortada en trozos de **≈25 g**.
- Determinación de carga inicial mediante diluciones sucesivas y siembra en superficie (método de extensión), tomando 1 mL de la serie 10⁻⁵ sobre placa con medio sólido específico por microorganismo.

**Resultado/conclusión que declara la tesis:**

- La carga microbiana inicial fue del orden de **2.3–2.4 × 10⁷ UFC/g**.
- Tras el tratamiento UV-C, la carga se redujo al orden de **0.6–0.7 × 10⁴ UFC/g**, logrando ≈ **97 %** de inactivación.
- La **mejor condición** de tratamiento es **3 minutos de exposición a 10 cm** de distancia entre la fuente UV-C y la muestra.
- La tesis clasifica su diseño como **preexperimental**.

---

## 3. Estructura experimental y factores

Aunque la tesis se declara **preexperimental** (esquema antes–después `O₁ X O₂`), la estructura real de los datos es **factorial con repeticiones**:

| Factor | Niveles | Tipo de factor |
|---|---|---|
| **Microorganismo** | *E. coli*, *S. typhimurium* | Factor biológico (2 niveles) |
| **Tiempo de exposición** | 1, 3 y 5 minutos | Factor experimental (3 niveles) |
| **Altura / distancia** | 5, 10 y 20 cm | Factor experimental (3 niveles) |
| **Repetición** | M1, M2, M3, M4 | Error experimental (4 réplicas) |

- **Variable respuesta (Y):** recuento de microorganismos sobrevivientes en placa Petri tras el tratamiento UV-C (UFC).
- Para datos microbiológicos se recomienda además trabajar en escala logarítmica:
  - `Y = log₁₀(UFC/g)`, o
  - **reducción logarítmica** `R = log₁₀(N₀) − log₁₀(N)`, con `N₀` = carga inicial y `N` = carga post-tratamiento.

**Total de observaciones:** 2 microorganismos × 3 tiempos × 3 alturas × 4 repeticiones = **72 datos** (más los controles `T(0)–H(0)`).

---

## 4. Datos crudos (recuentos en placa)

> Cada celda es el recuento de sobrevivientes; las 4 filas por tratamiento son las repeticiones M1–M4. La fila *Promedio* es la media de las 4 réplicas.

### 4.1. *Escherichia coli*

**Control** `T(0)–H(0)`: 225, 240, 245, 230 → **media = 235**

| Tiempo \ Altura | H1 = 5 cm | H2 = 10 cm | H3 = 20 cm |
|---|---|---|---|
| **T1 = 1 min** | 10, 8, 9, 11 | 12, 13, 12, 11 | 14, 17, 14, 15 |
| **T2 = 3 min** | 12, 13, 11, 11 | **6, 8, 7, 6** | 14, 12, 13, 15 |
| **T3 = 5 min** | 17, 18, 16, 15 | 15, 12, 13, 14 | 11, 13, 12, 14 |

| Promedios | H1 = 5 cm | H2 = 10 cm | H3 = 20 cm |
|---|---|---|---|
| **1 min** | 9.50 | 12.00 | 15.00 |
| **3 min** | 11.75 | **6.75** | 13.50 |
| **5 min** | 16.50 | 13.50 | 12.50 |

### 4.2. *Salmonella typhimurium*

**Control** `T(0)–H(0)`: 247, 245, 236, 233 → **media = 240.25**

| Tiempo \ Altura | H1 = 5 cm | H2 = 10 cm | H3 = 20 cm |
|---|---|---|---|
| **T1 = 1 min** | 13, 10, 12, 11 | 16, 17, 15, 17 | 17, 19, 20, 18 |
| **T2 = 3 min** | 10, 11, 10, 14 | **6, 5, 7, 6** | 15, 15, 16, 13 |
| **T3 = 5 min** | 17, 17, 19, 15 | 11, 15, 12, 14 | 11, 13, 10, 14 |

| Promedios | H1 = 5 cm | H2 = 10 cm | H3 = 20 cm |
|---|---|---|---|
| **1 min** | 11.50 | 16.25 | 18.50 |
| **3 min** | 11.25 | **6.00** | 14.75 |
| **5 min** | 17.00 | 13.00 | 12.00 |

> En **negrita** la combinación **3 min – 10 cm**, que arroja el menor recuento promedio en ambos microorganismos (6.75 y 6.00), coincidiendo con la "mejor condición" declarada por la tesis.

---

## 5. Inconsistencias detectadas en el libro de datos original

Durante la lectura del archivo `radiación UV-C para inactivar.xlsx` se detectaron incoherencias que el examen pide corregir (recomendación de "corregir inconsistencias entre tablas de recuento, promedios y porcentajes de inactivación"):

1. **Controles intercambiados.** En las hojas de recuento, el control de *E. coli* promedia **235** y el de *Salmonella* **240.25**. En la hoja resumen *Resultados de Inactivación* aparecen invertidos (E. coli = 240, Salmonella = 235), y el resumen de UFC/g del abstract (E. coli 2.3×10⁷ vs *Salmonella* 2.4×10⁷) tampoco coincide con esa hoja.
2. **Columna "% Reducción" de *Salmonella* mal calculada** en la hoja `Hoja8`: arroja valores del orden de 0.04–0.07 (no normalizados a 100) en lugar de ≈ 93–97 %, por un error en la fórmula de la celda.
3. La carga inicial reportada como UFC/g (≈10⁷) y el recuento en placa (≈10²) corresponden a escalas distintas (dilución 10⁻⁵), lo que debe explicitarse al interpretar el % de inactivación.

Estas observaciones **no afectan** el análisis de los 72 recuentos experimentales (que son consistentes), pero sí la presentación de resultados resumidos de la tesis.

---

## 6. Archivos fuente

| Archivo | Contenido |
|---|---|
| `enunciado/Presentacion  4 examen.pptx` | Consigna del examen (preguntas del docente). |
| `enunciado/radiación UV-C para inactivar...docx` | Tesis original (resumen, introducción, metodología). |
| `enunciado/radiación UV-C para inactivar.xlsx` | Datos crudos: hojas *Recuento E. coli*, *Recuento Salmonella*, *Hoja8*, *Resultados de Inactivacion*. |

---

## 7. Entregables de la resolución

| Entregable | Ruta |
|---|---|
| Documentación del enunciado | `enunciado.md` (este archivo) |
| Datos ordenados (tidy) | `codigo/datos_uvc.csv` |
| Notebook de análisis reproducible | `codigo/analisis_uvc.ipynb` |
| Figuras generadas | `codigo/figuras/*.png` |
| Informe estadístico final (LaTeX) | `resolucion/informe_mejorado.tex` (se compila a PDF) |
