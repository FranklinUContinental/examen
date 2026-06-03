# Radiación UV-C para la inactivación de *E. coli* y *S. typhimurium* en carne de cuy

> **Examen de Primera Unidad** · *Métodos Estadísticos para Investigaciones Experimentales*
> Escuela de Posgrado — Unidad de Posgrado en Ciencias Físicas y Matemáticas
> **Universidad Nacional de Trujillo** · **Docente:** Dr. Roger Reyna Segura
> **Estudiante:** Rudy Franklin Condori Quilla

![Status](https://img.shields.io/badge/estado-completo-success)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![LaTeX](https://img.shields.io/badge/documento-LaTeX%20%2F%20PDF-008080)
![Reproducible](https://img.shields.io/badge/an%C3%A1lisis-reproducible-brightgreen)

Este repositorio reúne la **auditoría estadística** de la tesis *«Radiación UV-C para inactivar
patógenos alimentarios* Escherichia coli *y* Salmonella typhimurium *en carne de cuy contaminado»*.
El objetivo es **revisar su diseño declarado**, **proponer el modelo
adecuado**, **validar supuestos**, **rehacer el ANOVA** y **documentar** el resultado final de
forma profesional en **Python**, **Markdown** y **LaTeX/PDF**.

### ¿Qué encontrarás aquí?

- `codigo/analisis_uvc.ipynb`: notebook reproducible con los modelos, tablas y figuras.
- `codigo/build_notebook.py`: script que reconstruye el notebook celda por celda.
- `resolucion/informe_mejorado.tex`: informe final en LaTeX (con carátula), listo para compilar a PDF.
- `codigo/figuras/`: figuras generadas por el análisis.

---

## 📌 Resumen ejecutivo

La siguiente síntesis responde, de forma directa, a las preguntas centrales del examen:

| Pregunta del examen | Respuesta de la auditoría |
|---|---|
| ¿El diseño y los resultados son correctos? | El **resultado experimental es aceptable** (~95 % de inactivación), pero el **diseño declarado «preexperimental» es inadecuado**. |
| ¿Cuál es el mejor diseño? | Un **DCA con arreglo factorial 2×3×3** (microorganismo × tiempo × altura), con 4 repeticiones (N = 72). |
| ¿Se cumplen los supuestos? | Homocedasticidad holgada (Levene); normalidad global aceptable (Shapiro–Wilk), dudosa en *E. coli* → se recomienda transformación logarítmica. |
| ¿Qué dice el nuevo análisis? | Son significativos microorganismo, tiempo y altura; **domina la interacción tiempo×altura** (F = 52.76). |
| ¿Conclusión final? | El **tratamiento óptimo es 3 min a 10 cm** (recuentos 6.75 y 6.00 UFC), confirmado por Tukey. |

---

## 📂 Estructura del proyecto

```
examen/
├── README.md                     ← este documento (guía rápida + resumen)
├── enunciado.md                  ← consigna del examen y datos crudos documentados
│
├── enunciado/                    ← material fuente original (sin modificar)
│   ├── Presentacion 4 examen.pptx                 · consigna del docente
│   ├── radiación UV-C ... contaminado.docx        · tesis original
│   └── radiación UV-C para inactivar.xlsx         · libro de datos crudos
│
├── codigo/                       ← análisis reproducible
│   ├── analisis_uvc.ipynb        · notebook con modelos, supuestos y figuras inline
│   ├── build_notebook.py         · script que reconstruye el notebook (nbformat)
│   ├── datos_uvc.csv             · datos ordenados (formato tidy, 72 obs.)
│   ├── requirements.txt          · dependencias de Python
│   └── figuras/                  · gráficos generados (01–06, PNG)
│
└── resolucion/                   ← informe final (documentación principal)
    └── informe_mejorado.tex      · ★ informe final en LaTeX con carátula (se compila a PDF)
```

---

## 🧭 Índice de la resolución

La resolución completa está documentada en LaTeX en
[`resolucion/EXAMEN.pdf`](resolucion/EXAMEN.pdf) . El
documento sigue, para cada modelo, la misma lógica de exposición: **pasos de trabajo**, **modelo
estadístico**, **script de Python**, **resultados** e **interpretación**. La estructura de
secciones del informe es la siguiente:

| § | Sección |
|---|---|
| 1 | Introducción y objetivos |
| 2 | Datos y estructura experimental |
| 3 | Evaluación crítica del diseño original |
| 4 | Estadística descriptiva y exploración gráfica |
| 5 | Modelo 1 — Diseño preexperimental |
| 6 | Modelo 2 — Diseño Completamente al Azar (DCA) |
| 7 | Modelo 3 — Diseño de Bloques Completamente al Azar (DBCA) |
| 8 | Modelo 4 — Cuadrado Latino (no aplica) |
| 9 | Modelo 5 — Factorial 3×3 por microorganismo |
| 10 | Modelo 6 — Factorial 2×3×3 (recomendado) |
| 11 | Validación de supuestos |
| 12 | Comparaciones múltiples (Tukey) |
| 13 | Comparación general de modelos |
| 14 | Discusión |
| 15 | Conclusiones y recomendaciones |
| 16 | Respuesta final tipo examen |

---

## 🔬 Resolución (síntesis)

Cada modelo del informe se presenta en **tres capas**: **(a)** resolución manual (modelo,
hipótesis y fórmulas), **(b)** script de Python que reproduce el cálculo y **(c)** resultados con
su interpretación. En el notebook y en el informe, cada bloque sigue además una secuencia clara:
**Paso 1** preparar o resumir los datos, **Paso 2** ejecutar el análisis o gráfico, y **Paso 3**
interpretar el resultado.

### 1. El diseño declarado es inadecuado
La tesis se declara **preexperimental** (esquema antes–después `O₁ X O₂`), pero la estructura real
de los datos es **factorial con repeticiones**. El preexperimental no permite estimar la
interacción tiempo×altura ni separar los efectos de los factores experimentales. **Veredicto:**
debió declararse un **DCA con arreglo factorial 2×3×3**.

### 2. Recorrido por los modelos
| Modelo | Aplicabilidad | Resultado | Evaluación |
|---|---|---|---|
| Preexperimental | Parcial | Reducción global ~95 % (t de Welch) | Insuficiente |
| DCA | Sí | Tratamientos significativos (p < 0.001) | Válido, no separa efectos |
| DBCA | Condicional | Bloques **no** significativos (p = 0.59/0.76) | No mejora al DCA |
| Cuadrado latino | **No** | Confundiría la interacción con el error | Inaplicable |
| Factorial 3×3 (por microorg.) | Sí | Tiempo, altura e interacción significativos | Muy adecuado |
| **Factorial 2×3×3** | **Sí** | Todos los efectos + interacciones clave | **★ Recomendado** |

### 3. ANOVA del modelo recomendado (factorial 2×3×3)
R² = 0.885; el término **dominante** es la interacción **tiempo×altura** (F = 52.76, p < 0.001).
La interacción microorganismo×altura y la triple **no** son significativas.

### 4. Supuestos y comparaciones múltiples
Homocedasticidad satisfecha (Levene p > 0.7); normalidad global aceptable (Shapiro–Wilk p = 0.068),
dudosa en *E. coli* (recomendable transformación log). **Tukey** confirma que **3 min – 10 cm** es
el tratamiento de menor recuento y difiere significativamente de la mayoría.

### 5. Conclusión
La radiación UV-C **sí** inactiva los patógenos (~95 %). La **conclusión experimental de la tesis
es aceptable**, pero su **análisis estadístico debe reformularse** como factorial 2×3×3. El
**tratamiento óptimo es 3 minutos a 10 cm** en ambos microorganismos.

> Desarrollo completo, ecuaciones, tablas de ANOVA y figuras: ver el informe en LaTeX
> [`resolucion/informe_mejorado.tex`](resolucion/informe_mejorado.tex) (compílalo para obtener el PDF).

---

## ▶️ Reproducir el análisis

Si quieres repetir el flujo de trabajo desde cero, sigue estos pasos:

```bash
# 1) Instalar dependencias
cd codigo
pip install -r requirements.txt

# 2) Ejecutar el notebook (genera resultados y figuras en codigo/figuras/)
jupyter notebook analisis_uvc.ipynb

# (opcional) Reconstruir el notebook desde el script fuente
python build_notebook.py
```

Todos los números del informe provienen de [`codigo/analisis_uvc.ipynb`](codigo/analisis_uvc.ipynb)
y han sido verificados contra los datos crudos del Excel original.

## 📄 Compilar la documentación LaTeX

Para obtener el PDF final, compila dos veces el archivo principal:

```bash
cd resolucion
pdflatex informe_mejorado.tex      # ejecutar dos veces para resolver referencias e índice
pdflatex informe_mejorado.tex
```

Requiere una distribución LaTeX (TeX Live / MiKTeX) con los paquetes `babel` (español),
`booktabs`, `graphicx`, `amsmath`, `hyperref` y `listings`. La ruta de figuras ya está
configurada con `\graphicspath{{../codigo/figuras/}{figuras/}}`.

---

## 📊 Figuras generadas

Las figuras se guardan en `codigo/figuras/` y se muestran también dentro del notebook.

| Figura | Descripción |
|---|---|
| `01_boxplots_factores.png` | Recuento por tiempo y por altura |
| `02_heatmaps_medias.png` | Mapa de calor de medias (mínimo en 3 min–10 cm) |
| `03_pct_inactivacion.png` | % de inactivación por tratamiento (todos > 94 %) |
| `04_interaccion_tiempo_altura.png` | Gráficos de interacción tiempo×altura |
| `05_diagnostico_residuos.png` | Q–Q plot y residuos vs. ajustados |
| `06_medias_tukey.png` | Medias por tratamiento con IC 95 % (Tukey) |

---

## 📚 Referencias

- Montgomery, D. C. (2019). *Design and Analysis of Experiments*. Wiley.
- Kuehl, R. O. (2001). *Diseño de experimentos*. Thomson.
- Steel, R. G. D., Torrie, J. H. & Dickey, D. A. (1997). *Principles and Procedures of Statistics*. McGraw-Hill.
- Gutiérrez Pulido, H. & De la Vara Salazar, R. (2012). *Análisis y diseño de experimentos*. McGraw-Hill.

---

<sub>Curso de *Métodos Estadísticos para Investigaciones Experimentales* — Escuela de Posgrado, Universidad Nacional de Trujillo.</sub>
