# -*- coding: utf-8 -*-
"""Construye codigo/analisis_uvc.ipynb celda por celda con nbformat."""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell
from pathlib import Path

cells = []
def md(s):   cells.append(new_markdown_cell(s))
def code(s): cells.append(new_code_cell(s))

# ----------------------------------------------------------------------
md(r"""# Análisis estadístico — Radiación UV-C en carne de cuy

**Reanálisis del experimento de inactivación de *E. coli* y *S. typhimurium* mediante radiación UV-C.**

Curso: *Métodos Estadísticos para Investigaciones Experimentales* — Escuela de Posgrado, UNT.

Este notebook reproduce de forma **íntegra y verificable** todos los modelos del informe:
preexperimental, DCA, DBCA, discusión del cuadrado latino, factorial 3×3 por microorganismo
y factorial 2×3×3, junto con la validación de supuestos (Shapiro–Wilk, Levene), las
comparaciones múltiples de Tukey y todos los gráficos.

Cada sección sigue el mismo esquema de trabajo: **(1)** preparar o resumir los datos,
**(2)** ejecutar el modelo o gráfico correspondiente y **(3)** interpretar el resultado.
Las figuras se guardan en disco y además se muestran dentro del notebook para facilitar su
revisión.

> Diseño recomendado: **DCA con arreglo factorial 2×3×3**. Variable respuesta: recuento de
> sobrevivientes en placa (UFC). Factores: microorganismo (2), tiempo (1, 3, 5 min) y altura
> (5, 10, 20 cm), con 4 repeticiones.""")

# ----------------------------------------------------------------------
md("## 1. Configuración del entorno")
code(r"""import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.graphics.factorplots import interaction_plot
import os, warnings
from IPython import get_ipython
from IPython.display import display, Image
from statsmodels.stats.multicomp import pairwise_tukeyhsd

warnings.filterwarnings("ignore")

sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams["figure.dpi"] = 110
FIG = "figuras_"
os.makedirs(FIG, exist_ok=True)
pd.set_option("display.float_format", lambda x: f"{x:,.4f}")
ip = get_ipython()
if ip is not None:
    try:
        ip.run_line_magic("matplotlib", "inline")
    except Exception:
        pass

def guardar_y_mostrar(fig, nombre):
    ruta = os.path.join(FIG, nombre)
    fig.savefig(ruta, bbox_inches="tight")
    plt.close(fig)
    display(Image(filename=ruta))
    return ruta

print("Entorno listo. Figuras ->", os.path.abspath(FIG))""")

# ----------------------------------------------------------------------
md("""## 2. Carga y ordenamiento de los datos (formato *tidy*)

**Paso 1.** Incrustar los recuentos en el notebook para que sea **autocontenido y reproducible**.

**Paso 2.** Organizar cada combinación tiempo×altura con 4 réplicas (M1–M4) en una tabla ordenada.

**Paso 3.** Guardar los controles `T(0)–H(0)` para el modelo preexperimental y para calcular la
reducción respecto a `N0`.""")
code(r"""# Recuentos de sobrevivientes en placa (UFC). Orden de réplicas: M1, M2, M3, M4
tratados = {
 "E. coli": {
   (1,5):[10,8,9,11],  (1,10):[12,13,12,11], (1,20):[14,17,14,15],
   (3,5):[12,13,11,11], (3,10):[6,8,7,6],     (3,20):[14,12,13,15],
   (5,5):[17,18,16,15], (5,10):[15,12,13,14], (5,20):[11,13,12,14]},
 "S. typhimurium": {
   (1,5):[13,10,12,11], (1,10):[16,17,15,17], (1,20):[17,19,20,18],
   (3,5):[10,11,10,14], (3,10):[6,5,7,6],     (3,20):[15,15,16,13],
   (5,5):[17,17,19,15], (5,10):[11,15,12,14], (5,20):[11,13,10,14]},
}
control = {"E. coli":[225,240,245,230], "S. typhimurium":[247,245,236,233]}

filas = []
for org, d in tratados.items():
    for (t,h), vals in d.items():
        for rep, v in enumerate(vals, 1):
            filas.append(dict(microorganismo=org, tiempo=t, altura=h, rep=rep, recuento=v))
df = pd.DataFrame(filas)
df["trat"] = "T"+df.tiempo.astype(str)+"_H"+df.altura.astype(str)
df.to_csv("datos_uvc.csv", index=False)
print("Observaciones tratadas:", len(df), "| controles:", sum(len(v) for v in control.values()))
df.head(10)""")

# ----------------------------------------------------------------------
md("""## 3. Estadística descriptiva

**Paso 1.** Calcular `N0` como la media de los controles por microorganismo.

**Paso 2.** Resumir cada tratamiento con media, desviación estándar y número de réplicas.

**Paso 3.** Estimar la reducción logarítmica `R = log10(N0) − log10(N)` y el porcentaje de
inactivación para interpretar la magnitud del efecto.""")
code(r"""N0 = {org: np.mean(v) for org, v in control.items()}
print("Carga media del control (N0):", {k: round(v,2) for k,v in N0.items()})

resumen = (df.groupby(["microorganismo","tiempo","altura"])["recuento"]
             .agg(media="mean", sd="std", n="count").reset_index())
resumen["reduccion_log10"] = resumen.apply(
    lambda r: np.log10(N0[r.microorganismo]) - np.log10(r.media), axis=1)
resumen["pct_inactivacion"] = resumen.apply(
    lambda r: 100*(1 - r.media/N0[r.microorganismo]), axis=1)
resumen.round(3)""")

code(r"""# Tabla de medias en formato cruzado (tiempo x altura) por microorganismo
for org in df.microorganismo.unique():
    piv = (df[df.microorganismo==org]
           .pivot_table(index="tiempo", columns="altura", values="recuento", aggfunc="mean"))
    print(f"\nMedias de recuento — {org}")
    print(piv.round(2))
    mejor = resumen[resumen.microorganismo==org].sort_values("media").iloc[0]
    print(f"  -> menor recuento: {mejor.media:.2f} en tiempo={int(mejor.tiempo)} min, altura={int(mejor.altura)} cm")""")

# ----------------------------------------------------------------------
md("""## 4. Exploración gráfica

**Paso 1.** Explorar la distribución de recuentos por tiempo y por altura.

**Paso 2.** Comparar visualmente ambas especies para anticipar diferencias entre factores.

### 4.1. Distribución de recuentos por factor""")
code(r"""fig, axes = plt.subplots(1, 2, figsize=(12,4.5))
sns.boxplot(data=df, x="tiempo", y="recuento", hue="microorganismo", ax=axes[0])
axes[0].set(title="Recuento por tiempo de exposición", xlabel="Tiempo (min)", ylabel="Recuento (UFC)")
sns.boxplot(data=df, x="altura", y="recuento", hue="microorganismo", ax=axes[1])
axes[1].set(title="Recuento por altura de irradiación", xlabel="Altura (cm)", ylabel="Recuento (UFC)")
plt.tight_layout(); guardar_y_mostrar(fig, "01_boxplots_factores.png")""")

md("### 4.2. Mapas de calor de medias (tiempo × altura)")
code(r"""fig, axes = plt.subplots(1, 2, figsize=(11,4.2))
for ax, org in zip(axes, df.microorganismo.unique()):
    piv = df[df.microorganismo==org].pivot_table(index="tiempo", columns="altura",
                                                  values="recuento", aggfunc="mean")
    sns.heatmap(piv, annot=True, fmt=".2f", cmap="RdYlGn_r", ax=ax, cbar_kws={"label":"Recuento medio"})
    ax.set(title=f"{org}", xlabel="Altura (cm)", ylabel="Tiempo (min)")
plt.suptitle("Recuento medio de sobrevivientes (menor = mejor inactivación)", y=1.03)
plt.tight_layout(); guardar_y_mostrar(fig, "02_heatmaps_medias.png")""")

md("""### 4.3. % de inactivación por tratamiento

**Paso 3.** Resumir el porcentaje de inactivación por combinación tiempo×altura y validar que todas
las condiciones superan el 94 %.""")
code(r"""orden = resumen.assign(lbl=resumen.tiempo.astype(str)+"min/"+resumen.altura.astype(str)+"cm")
fig, ax = plt.subplots(figsize=(11,4.5))
sns.barplot(data=orden, x="lbl", y="pct_inactivacion", hue="microorganismo", ax=ax)
ax.set(title="% de inactivación respecto al control", xlabel="Tiempo / Altura", ylabel="% inactivación")
ax.set_ylim(90, 100); ax.tick_params(axis="x", rotation=30)
plt.tight_layout(); guardar_y_mostrar(fig, "03_pct_inactivacion.png")""")

# ----------------------------------------------------------------------
md(r"""## 5. Modelo 1 — Diseño preexperimental

**Paso 1.** Comparar la condición **control** (sin tratar) contra la condición **tratada**
(UV-C aplicado), es decir el esquema antes–después \(O_1\; X\; O_2\).

**Paso 2.** Formular la hipótesis nula y alternativa para ambos microorganismos.

**Paso 3.** Aplicar la prueba *t* de Welch para verificar si la radiación UV-C reduce la carga
microbiana.

$$Y_{ij} = \mu + \tau_i + e_{ij}$$

**Hipótesis:** \(H_0:\ \mu_{control} = \mu_{tratado}\) vs \(H_1:\ \mu_{control} \neq \mu_{tratado}\).

Es un *t* de dos muestras (control vs. tratado) por microorganismo. Sólo responde *si hubo
reducción*, no *cuál tiempo/altura es óptimo*.""")
code(r"""print("MODELO PREEXPERIMENTAL (control vs tratado)\n" + "="*48)
for org in df.microorganismo.unique():
    c = np.array(control[org], dtype=float)
    t = df[df.microorganismo==org]["recuento"].values.astype(float)
    tstat, p = stats.ttest_ind(c, t, equal_var=False)
    print(f"\n{org}:")
    print(f"  media control = {c.mean():.2f} | media tratado = {t.mean():.2f}")
    print(f"  reducción     = {100*(1-t.mean()/c.mean()):.2f}%")
    print(f"  t = {tstat:.3f}, p = {p:.3e}  -> {'RECHAZA H0 (hay efecto)' if p<0.05 else 'no rechaza H0'}")""")

# ----------------------------------------------------------------------
md(r"""## 6. Modelo 2 — Diseño Completamente al Azar (DCA)

**Paso 1.** Tratar cada combinación tiempo×altura como **un** tratamiento independiente
(9 tratamientos).

**Paso 2.** Ajustar un ANOVA de un factor por microorganismo.

**Paso 3.** Interpretar si existen diferencias significativas entre tratamientos, aunque sin
separar todavía tiempo, altura e interacción.

$$Y_{ij} = \mu + \tau_i + e_{ij}, \qquad H_0:\ \mu_1=\cdots=\mu_9$$

Se ajusta por microorganismo.""")
code(r"""def anova_dca(sub, org):
    modelo = ols("recuento ~ C(trat)", data=sub).fit()
    tabla = anova_lm(modelo, typ=2)
    print(f"\nANOVA DCA — {org}")
    print(tabla.round(4))
    return modelo

mod_dca = {}
for org in df.microorganismo.unique():
    mod_dca[org] = anova_dca(df[df.microorganismo==org], org)""")

# ----------------------------------------------------------------------
md(r"""## 7. Modelo 3 — Diseño de Bloques Completamente al Azar (DBCA)

**Paso 1.** Considerar la repetición (M1–M4) como **bloque** y verificar si aporta variación
sistemática.

**Paso 2.** Ajustar el modelo con tratamientos + bloques.

**Paso 3.** Evaluar el p-valor del bloque para decidir si el bloqueo es útil.

$$Y_{ij} = \mu + \tau_i + \beta_j + e_{ij}$$

Permite verificar si el bloqueo (réplica) aporta. Se espera que **no** sea significativo,
porque las réplicas no representan una fuente sistemática real de variación.""")
code(r"""def anova_dbca(sub, org):
    modelo = ols("recuento ~ C(trat) + C(rep)", data=sub).fit()
    tabla = anova_lm(modelo, typ=2)
    print(f"\nANOVA DBCA — {org}  (tratamientos + bloque=réplica)")
    print(tabla.round(4))
    p_bloque = tabla.loc["C(rep)","PR(>F)"]
    print(f"  Bloques: p = {p_bloque:.3f} -> {'significativo' if p_bloque<0.05 else 'NO significativo (bloqueo no aporta)'}")
    return modelo

mod_dbca = {}
for org in df.microorganismo.unique():
    mod_dbca[org] = anova_dbca(df[df.microorganismo==org], org)""")

# ----------------------------------------------------------------------
md(r"""## 8. Modelo 4 — Cuadrado Latino (¿aplica?)

**Paso 1.** Revisar si el cuadrado latino es compatible con la estructura del experimento.

**Paso 2.** Verificar si tiempo y altura actúan como bloques o como factores de interés.

**Paso 3.** Comparar el modelo aditivo con el modelo que incluye interacción para mostrar por
qué el cuadrado latino **no corresponde** a este caso.

El cuadrado latino exige **dos** factores de bloqueo cruzados y un número igual de
tratamientos, filas y columnas. **No corresponde** a este experimento porque:

- Tiempo y altura **no** son bloques: son **factores experimentales de interés**.
- Su **interacción** tiempo×altura es justamente lo que se quiere estimar, y el cuadrado
  latino la confunde con el error.
- No hay dos fuentes de ruido cruzadas que controlar.

A continuación se demuestra **numéricamente** por qué descartarlo: si se forzara un modelo
puramente aditivo (sin interacción), la interacción quedaría escondida en el error.""")
code(r"""org = "E. coli"
sub = df[df.microorganismo==org]
aditivo = ols("recuento ~ C(tiempo) + C(altura)", data=sub).fit()
completo = ols("recuento ~ C(tiempo) * C(altura)", data=sub).fit()
print(f"{org}: modelo ADITIVO (estilo cuadrado latino, sin interacción)")
print(f"  SC error aditivo  = {aditivo.ssr:.2f}  (gl={int(aditivo.df_resid)})")
print(f"  SC error completo = {completo.ssr:.2f}  (gl={int(completo.df_resid)})")
print(f"  -> La interacción tiempo×altura aporta SC = {aditivo.ssr-completo.ssr:.2f},")
print("     que un diseño tipo cuadrado latino enviaría erróneamente al error.")""")

# ----------------------------------------------------------------------
md(r"""## 9. Modelo 5 — Factorial 3×3 por microorganismo

**Paso 1.** Separar tiempo y altura como factores experimentales de interés.

**Paso 2.** Ajustar el factorial 3×3 por microorganismo para descomponer el efecto total.

**Paso 3.** Interpretar los efectos principales y la interacción tiempo×altura.

$$Y_{ijk} = \mu + \alpha_i + \beta_j + (\alpha\beta)_{ij} + e_{ijk}$$

con \(\alpha\) = tiempo, \(\beta\) = altura. Descompone los 9 tratamientos del DCA en
**efectos principales** (tiempo, altura) e **interacción**.""")
code(r"""mod_fac = {}
for org in df.microorganismo.unique():
    sub = df[df.microorganismo==org]
    modelo = ols("recuento ~ C(tiempo) * C(altura)", data=sub).fit()
    tabla = anova_lm(modelo, typ=2)
    mod_fac[org] = modelo
    print(f"\nANOVA Factorial 3x3 — {org}")
    print(tabla.round(4))
    print(f"  R^2 = {modelo.rsquared:.3f}")""")

md("### 9.1. Gráficos de interacción tiempo × altura")
code(r"""fig, axes = plt.subplots(1, 2, figsize=(12,4.5))
for ax, org in zip(axes, df.microorganismo.unique()):
    sub = df[df.microorganismo==org]
    interaction_plot(x=sub.tiempo.values, trace=sub.altura.values, response=sub.recuento.values,
                     ax=ax, markers=["o","s","^"], colors=["#2c7fb8","#d95f02","#1b9e77"])
    ax.set(title=f"Interacción tiempo×altura — {org}", xlabel="Tiempo (min)", ylabel="Recuento medio")
    ax.legend(title="Altura (cm)")
plt.tight_layout(); guardar_y_mostrar(fig, "04_interaccion_tiempo_altura.png")""")

# ----------------------------------------------------------------------
md(r"""## 10. Modelo 6 — Factorial 2×3×3 (modelo recomendado)

**Paso 1.** Integrar microorganismo, tiempo y altura en un único modelo factorial.

**Paso 2.** Evaluar todas las interacciones posibles para no perder información.

**Paso 3.** Usar este modelo como referencia principal porque representa la estructura real del
experimento.

$$Y_{ijkl} = \mu + A_i + B_j + C_k + (AB)_{ij} + (AC)_{ik} + (BC)_{jk} + (ABC)_{ijk} + e_{ijkl}$$

con \(A\) = microorganismo, \(B\) = tiempo, \(C\) = altura. Es el modelo **más completo**:
compara microorganismos y todas las interacciones en un solo análisis.""")
code(r"""modelo_full = ols("recuento ~ C(microorganismo) * C(tiempo) * C(altura)", data=df).fit()
tabla_full = anova_lm(modelo_full, typ=2)
print("ANOVA Factorial 2x3x3 (modelo completo)")
print(tabla_full.round(4))
print(f"\nR^2 = {modelo_full.rsquared:.3f} | R^2 ajustado = {modelo_full.rsquared_adj:.3f}")
tabla_full""")

# ----------------------------------------------------------------------
md(r"""## 11. Validación de supuestos del modelo 2×3×3

**Paso 1.** Extraer los residuos del modelo completo.

**Paso 2.** Evaluar normalidad con Shapiro–Wilk y homogeneidad con Levene.

**Paso 3.** Complementar con diagnóstico gráfico de residuos.

Sobre los residuos del modelo completo:
- **Normalidad:** Shapiro–Wilk (global y por microorganismo).
- **Homogeneidad de varianzas:** Levene entre los 9 tratamientos.
- **Diagnóstico gráfico:** Q–Q plot y residuos vs. ajustados.""")
code(r"""res = modelo_full.resid
print("NORMALIDAD (Shapiro-Wilk)")
W, p = stats.shapiro(res)
print(f"  Global: W = {W:.4f}, p = {p:.4f} -> {'normal' if p>0.05 else 'dudosa'}")
for org in df.microorganismo.unique():
    m = ols("recuento ~ C(tiempo)*C(altura)", data=df[df.microorganismo==org]).fit()
    W, p = stats.shapiro(m.resid)
    print(f"  {org}: W = {W:.4f}, p = {p:.4f} -> {'normal' if p>0.05 else 'dudosa'}")

print("\nHOMOGENEIDAD DE VARIANZAS (Levene, 9 tratamientos)")
grupos = [g["recuento"].values for _,g in df.groupby(["microorganismo","trat"])]
W, p = stats.levene(*grupos, center="median")
print(f"  Global: estadístico = {W:.4f}, p = {p:.4f} -> {'homogéneas' if p>0.05 else 'heterogéneas'}")
for org in df.microorganismo.unique():
    gs = [g["recuento"].values for _,g in df[df.microorganismo==org].groupby("trat")]
    W, p = stats.levene(*gs, center="median")
    print(f"  {org}: estadístico = {W:.4f}, p = {p:.4f} -> {'homogéneas' if p>0.05 else 'heterogéneas'}")""")

code(r"""fig, axes = plt.subplots(1, 2, figsize=(12,4.5))
sm.qqplot(res, line="s", ax=axes[0])
axes[0].set_title("Q–Q plot de residuos (modelo 2×3×3)")
axes[1].scatter(modelo_full.fittedvalues, res, alpha=0.7, edgecolor="k", linewidth=0.3)
axes[1].axhline(0, color="red", ls="--")
axes[1].set(title="Residuos vs. valores ajustados", xlabel="Ajustado", ylabel="Residuo")
plt.tight_layout(); guardar_y_mostrar(fig, "05_diagnostico_residuos.png")""")

# ----------------------------------------------------------------------
md(r"""## 12. Comparaciones múltiples (Tukey HSD)

**Paso 1.** Confirmar que la interacción tiempo×altura es significativa.

**Paso 2.** Aplicar Tukey sobre los 9 tratamientos por microorganismo.

**Paso 3.** Identificar el tratamiento óptimo por comparación múltiple de medias.

Como la interacción tiempo×altura es significativa, el tratamiento óptimo se identifica por
**combinación de niveles**. Se aplica Tukey sobre los 9 tratamientos por microorganismo.""")
code(r"""for org in df.microorganismo.unique():
    sub = df[df.microorganismo==org]
    tuk = pairwise_tukeyhsd(sub.recuento, sub.trat, alpha=0.05)
    print(f"\nTukey HSD — {org}")
    print(tuk.summary())
    medias = sub.groupby("trat").recuento.mean().sort_values()
    print(f"  Tratamiento con menor recuento: {medias.index[0]} (media={medias.iloc[0]:.2f})")""")

code(r"""# Visualización de las medias por tratamiento con IC
fig, axes = plt.subplots(1, 2, figsize=(13,4.5), sharey=True)
for ax, org in zip(axes, df.microorganismo.unique()):
    sub = df[df.microorganismo==org]
    orden = sub.groupby("trat").recuento.mean().sort_values().index
    sns.pointplot(data=sub, x="trat", y="recuento", order=orden, ax=ax,
                  errorbar="ci", capsize=.3, color="#2c7fb8")
    ax.set(title=f"Medias por tratamiento (IC95%) — {org}", xlabel="Tratamiento (Ttiempo_Haltura)",
           ylabel="Recuento")
    ax.tick_params(axis="x", rotation=45)
plt.tight_layout(); guardar_y_mostrar(fig, "06_medias_tukey.png")""")

# ----------------------------------------------------------------------
md("""## 13. Comparación general de modelos

**Paso 1.** Reunir, en una sola tabla, la utilidad de cada modelo evaluado.

**Paso 2.** Comparar su aplicabilidad y su capacidad para explicar la variación observada.

**Paso 3.** Seleccionar el modelo recomendado a partir de la evidencia previa.""")
code(r"""comparacion = pd.DataFrame([
 ["Preexperimental","Parcial","Detecta reducción global","Insuficiente: no evalúa tiempo/altura/interacción"],
 ["DCA","Sí","Tratamientos significativos (p<0.001)","Válido, pero no separa efectos principales"],
 ["DBCA","Condicional","Tratamientos signif.; bloques NO","No mejora al DCA; bloqueo injustificado"],
 ["Cuadrado latino","No","No corresponde a la estructura","Confundiría la interacción con el error"],
 ["Factorial 3×3","Sí","Tiempo, altura e interacción signif.","Muy adecuado por microorganismo"],
 ["Factorial 2×3×3","Sí","Evalúa microorg., tiempo, altura e interacc.","MODELO RECOMENDADO (más completo)"],
], columns=["Modelo","Aplicabilidad","Resultado principal","Evaluación"])
comparacion""")

# ----------------------------------------------------------------------
md(r"""## 14. Conclusiones del reanálisis

**Paso 1.** Integrar la evidencia obtenida en los modelos y en los gráficos.

**Paso 2.** Sintetizar qué diseño es correcto y cuál es el tratamiento óptimo.

**Paso 3.** Redactar conclusiones y recomendaciones finales de forma clara y verificable.

1. **Hay efecto UV-C:** el modelo preexperimental confirma una reducción ≈95 % respecto al
   control en ambos patógenos (*p* ≪ 0.001).
2. **El diseño correcto es factorial**, no preexperimental: el experimento cruza tiempo y
   altura con réplicas.
3. **DCA y DBCA** detectan diferencias entre tratamientos, pero el bloqueo por réplica **no
   es significativo** → el DBCA no aporta sobre el DCA.
4. **El cuadrado latino no aplica:** tiempo y altura son factores, no bloques.
5. **El factorial 2×3×3 es el modelo recomendado.** Resultan significativos microorganismo,
   tiempo, altura, la interacción microorganismo×tiempo y, sobre todo, la **interacción
   tiempo×altura** (*p* ≪ 0.001). La interacción triple no es significativa.
6. **Supuestos:** homogeneidad de varianzas se cumple holgadamente; la normalidad global es
   aceptable (en *E. coli* es dudosa → se sugiere transformación log).
7. **Tratamiento óptimo (confirmado por Tukey):** **3 min – 10 cm**, con recuentos medios de
   **6.75** (*E. coli*) y **6.00** (*S. typhimurium*) — los más bajos del experimento.
8. **Veredicto:** la conclusión experimental de la tesis es correcta, pero su **diseño y
   análisis estadístico deben reformularse** (factorial, hipótesis bien planteadas, supuestos
   validados, ANOVA completo, Tukey y, preferiblemente, escala logarítmica).""")

# ----------------------------------------------------------------------
nb = new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name":"Python 3","language":"python","name":"python3"}
nb.metadata["language_info"] = {"name":"python"}
BASE_DIR = Path(__file__).resolve().parent
with open(BASE_DIR / "analisis_uvc.ipynb","w",encoding="utf-8") as f:
    nbf.write(nb, f)
print("Notebook creado:", BASE_DIR / "analisis_uvc.ipynb", "con", len(cells), "celdas")
