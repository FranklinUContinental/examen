# Resolución mediante diseños factoriales — efectos fijos, aleatorios y mixtos

> **Segunda resolución del examen** · *Métodos Estadísticos para Investigaciones Experimentales*
> Escuela de Posgrado — Unidad de Posgrado en Ciencias Físicas y Matemáticas
> **Universidad Nacional de Trujillo** · **Docente:** Dr. Roger Reyna Segura
> **Estudiante:** Rudy Franklin Condori Quilla

Esta carpeta resuelve el experimento de **radiación UV-C** (inactivación de *E. coli* y
*S. typhimurium* en carne de cuy) aplicando el **arreglo factorial 2×3×3** bajo las **tres
condiciones** del archivo [`condiciones.txt`](condiciones.txt):

1. **Modelo I — Efectos fijos**
2. **Modelo II — Efectos aleatorios**
3. **Modelo III — Efectos mixtos**

Sigue los mismos procedimientos de la primera resolución (carpetas `codigo/` y `resolucion/`,
estructura de tres capas por modelo: resolución manual → script de Python → resultados), selecciona
el **mejor modelo factorial** y lo **compara con la mejor opción de DCA** obtenida en el primer examen.

---

## 📌 Resumen ejecutivo

| Pregunta | Respuesta |
|---|---|
| ¿En qué se diferencian los tres modelos? | Comparten la **descomposición de sumas de cuadrados**; cambian los **cuadrados medios esperados (EMS)** y, por tanto, el **denominador de cada F**. |
| ¿Efectos fijos? | Significativos microorganismo, tiempo, altura, A×B y, dominante, la **interacción tiempo×altura** (F = 52.76). |
| ¿Efectos aleatorios? | Los efectos principales pierden significancia; aparecen **componentes de varianza negativos** → supuesto no coherente. |
| ¿Efectos mixtos (A fijo; B,C aleatorios)? | La **interacción tiempo×altura sigue muy significativa**; los efectos principales no, frente a CM_BC. |
| ¿Mejor modelo factorial? | **Efectos fijos**: niveles elegidos deliberadamente; máxima potencia (gl error = 54). |
| ¿Factorial vs. DCA? | El factorial fijo **contiene** al DCA (misma SSE y R²) y además descompone la variación → **superior**. |
| ¿Tratamiento óptimo? | **3 min – 10 cm** en ambos patógenos (recuentos 6.75 y 6.00). |

> **Hallazgo clave:** la interacción **tiempo×altura** es el único efecto significativo en las **tres**
> condiciones. Las conclusiones sobre los efectos principales **dependen del supuesto** fijos/aleatorios/mixtos.

---

## 📂 Estructura

```
resolucion_factorial/
├── README.md                    ← este documento
├── condiciones.txt              ← las 3 condiciones del modelo
│
├── codigo/
│   ├── analisis_factorial.ipynb · notebook con los 3 modelos, EMS, componentes de varianza y figuras
│   ├── analisis.py              · núcleo de cálculo (funciones puras y testeables)
│   ├── test_analisis.py         · pruebas unitarias (pytest) del núcleo de cálculo
│   ├── build_notebook.py        · script que reconstruye el notebook (nbformat)
│   ├── datos_uvc.csv            · datos ordenados (tidy, 72 obs.)
│   ├── requirements.txt         · dependencias de Python
│   └── figuras/                 · gráficos generados (PNG)
│
└── resolucion/
    └── informe_factorial.tex    · ★ informe final en LaTeX (se compila a PDF)
```

---

## 🧭 Índice del informe

| § | Sección |
|---|---|
| 1.1 | **Tabla resumen de los 6 modelos del primer examen** |
| 1–3 | Introducción, objetivos y estructura factorial |
| 4 | Marco teórico: EMS de efectos fijos, aleatorios y mixtos |
| 5 | ANOVA factorial base (descomposición común) |
| 6 | Modelo I — Efectos fijos |
| 7 | Modelo II — Efectos aleatorios (+ componentes de varianza) |
| 8 | Modelo III — Efectos mixtos |
| 9 | Comparación de las tres condiciones |
| 10 | Componentes de varianza e interacción |
| 11 | Selección del mejor modelo factorial |
| 12 | Comparación factorial vs. DCA del primer examen |
| 13–14 | Supuestos y Tukey |
| 15–16 | Conclusiones y respuesta final |

---

## ▶️ Reproducir el análisis

```bash
cd resolucion_factorial/codigo
pip install -r requirements.txt
python build_notebook.py                                    # (re)genera el notebook
jupyter nbconvert --to notebook --execute --inplace analisis_factorial.ipynb
```

## ✅ Pruebas unitarias

El núcleo de cálculo está aislado en `codigo/analisis.py` (funciones puras) y verificado con
`codigo/test_analisis.py` (**22 pruebas**, pytest). Validan: estructura del diseño, grados de
libertad y cuadrados medios, las pruebas F de los **tres modelos** (fijos, aleatorios con
Satterthwaite, y mixtos), los componentes de varianza, la equivalencia factorial = DCA descompuesto,
los supuestos y el tratamiento óptimo.

```bash
cd resolucion_factorial/codigo
pytest -v
```

## 📄 Compilar el informe LaTeX

```bash
cd resolucion_factorial/resolucion
pdflatex informe_factorial.tex     # ejecutar dos veces (índice y referencias)
pdflatex informe_factorial.tex
```

Requiere una distribución LaTeX (TeX Live / MiKTeX) con `babel` (español), `booktabs`, `graphicx`,
`amsmath`, `hyperref`, `longtable` y `listings`. La ruta de figuras ya apunta a `../codigo/figuras/`.

---

## 📊 Figuras

| Figura | Descripción |
|---|---|
| `03_F_por_condicion.png` | Valor F de cada efecto bajo fijos/aleatorios/mixtos (escala log) |
| `04_componentes_varianza.png` | Componentes de varianza (modelo aleatorio); domina tiempo×altura |
| `05_interaccion_tiempo_altura.png` | Gráficos de interacción tiempo×altura por microorganismo |
| `06_factorial_vs_dca.png` | Descomposición de la SC de tratamientos del DCA en el factorial |
| `07_diagnostico_residuos.png` | Q–Q plot y residuos vs. ajustados |

---

## 📚 Referencias

- Montgomery, D. C. (2019). *Design and Analysis of Experiments*. Wiley.
- Kuehl, R. O. (2001). *Diseño de experimentos*. Thomson.
- Steel, R. G. D., Torrie, J. H. & Dickey, D. A. (1997). *Principles and Procedures of Statistics*. McGraw-Hill.
- Gutiérrez Pulido, H. & De la Vara Salazar, R. (2012). *Análisis y diseño de experimentos*. McGraw-Hill.
