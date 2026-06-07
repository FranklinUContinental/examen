# -*- coding: utf-8 -*-
"""Núcleo de cálculo del análisis factorial UV-C bajo efectos fijos, aleatorios
y mixtos.

Reúne en funciones puras y reutilizables la misma lógica que desarrolla el
notebook ``analisis_factorial.ipynb``. Al estar aislada del notebook, puede
verificarse con pruebas unitarias (``test_analisis.py``).

Convenciones de los factores del arreglo factorial 2x3x3:
    A = microorganismo (a=2),  B = tiempo (b=3),  C = altura (c=3),  n=4 réplicas.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

# Constantes del diseño (con sufijo _NIV para no colisionar con la función C() de patsy)
A_NIV, B_NIV, C_NIV, N_REP = 2, 3, 3, 4
DATOS = Path(__file__).resolve().parent / "datos_uvc.csv"

# Nombre de cada término del modelo -> etiqueta de statsmodels
KEY = {
    "A": "C(microorganismo)",
    "B": "C(tiempo)",
    "C": "C(altura)",
    "AB": "C(microorganismo):C(tiempo)",
    "AC": "C(microorganismo):C(altura)",
    "BC": "C(tiempo):C(altura)",
    "ABC": "C(microorganismo):C(tiempo):C(altura)",
    "E": "Residual",
}


def cargar_datos(path: Path | str = DATOS) -> pd.DataFrame:
    """Lee los 72 recuentos en formato tidy."""
    return pd.read_csv(path)


def ajustar_factorial(df: pd.DataFrame):
    """Ajusta el modelo factorial completo 2x3x3 y devuelve (modelo, tabla_anova)."""
    modelo = ols("recuento ~ C(microorganismo)*C(tiempo)*C(altura)", data=df).fit()
    tabla = anova_lm(modelo, typ=2)  # balanceado: typ 1 = 2 = 3
    return modelo, tabla


def cuadrados_medios(tabla: pd.DataFrame):
    """Extrae SC, gl y CM por término a partir de la tabla ANOVA. Devuelve (SC, GL, MS)."""
    SC = {k: float(tabla.loc[v, "sum_sq"]) for k, v in KEY.items()}
    GL = {k: int(tabla.loc[v, "df"]) for k, v in KEY.items()}
    MS = {k: SC[k] / GL[k] for k in SC}
    return SC, GL, MS


# ---------------------------------------------------------------------------
# Pruebas F genéricas según el denominador que imponen los EMS
# ---------------------------------------------------------------------------
def f_simple(MS, GL, num: str, den: str):
    """F con denominador = un solo cuadrado medio. Devuelve (F, gl1, gl2, p)."""
    F = MS[num] / MS[den]
    p = float(stats.f.sf(F, GL[num], GL[den]))
    return F, GL[num], GL[den], p


def f_satterthwaite(MS, GL, num: str, mas, menos):
    """F con denominador sintético = sum(CM[mas]) - sum(CM[menos]).

    Devuelve (F, gl1, gl2_satterthwaite, p, denominador).
    """
    den = sum(MS[t] for t in mas) - sum(MS[t] for t in menos)
    gl_den = den ** 2 / sum((MS[t] ** 2) / GL[t] for t in (list(mas) + list(menos)))
    F = MS[num] / den
    p = float(stats.f.sf(F, GL[num], gl_den))
    return F, GL[num], gl_den, p, den


def modelo_fijos(MS, GL) -> dict:
    """Modelo I: todos los efectos contra el error."""
    return {k: f_simple(MS, GL, k, "E") for k in ["A", "B", "C", "AB", "AC", "BC", "ABC"]}


def modelo_aleatorios(MS, GL) -> dict:
    """Modelo II: denominadores según EMS del modelo de efectos aleatorios."""
    out = {
        "ABC": f_simple(MS, GL, "ABC", "E"),
        "AB": f_simple(MS, GL, "AB", "ABC"),
        "AC": f_simple(MS, GL, "AC", "ABC"),
        "BC": f_simple(MS, GL, "BC", "ABC"),
    }
    combos = {"A": (["AB", "AC"], ["ABC"]),
              "B": (["AB", "BC"], ["ABC"]),
              "C": (["AC", "BC"], ["ABC"])}
    for k, (mas, menos) in combos.items():
        F, g1, g2, p, _ = f_satterthwaite(MS, GL, k, mas, menos)
        out[k] = (F, g1, g2, p)
    return out


def modelo_mixtos(MS, GL) -> dict:
    """Modelo III: A (microorganismo) fijo; B (tiempo) y C (altura) aleatorios; restringido."""
    F_A, g1, g2, p, _ = f_satterthwaite(MS, GL, "A", ["AB", "AC"], ["ABC"])
    return {
        "A": (F_A, g1, g2, p),
        "B": f_simple(MS, GL, "B", "BC"),
        "C": f_simple(MS, GL, "C", "BC"),
        "AB": f_simple(MS, GL, "AB", "ABC"),
        "AC": f_simple(MS, GL, "AC", "ABC"),
        "BC": f_simple(MS, GL, "BC", "E"),
        "ABC": f_simple(MS, GL, "ABC", "E"),
    }


def componentes_varianza(MS, a=A_NIV, b=B_NIV, c=C_NIV, n=N_REP) -> dict:
    """Componentes de varianza por el método de los momentos (modelo aleatorio)."""
    return {
        "A": (MS["A"] - MS["AB"] - MS["AC"] + MS["ABC"]) / (b * c * n),
        "B": (MS["B"] - MS["AB"] - MS["BC"] + MS["ABC"]) / (a * c * n),
        "C": (MS["C"] - MS["AC"] - MS["BC"] + MS["ABC"]) / (a * b * n),
        "AB": (MS["AB"] - MS["ABC"]) / (c * n),
        "AC": (MS["AC"] - MS["ABC"]) / (b * n),
        "BC": (MS["BC"] - MS["ABC"]) / (a * n),
        "ABC": (MS["ABC"] - MS["E"]) / n,
        "error": MS["E"],
    }


def comparacion_dca(df: pd.DataFrame) -> dict:
    """Por microorganismo: DCA (9 tratamientos) vs. factorial 3x3.

    Devuelve, por organismo, las SC y el R^2, verificando que el factorial
    descompone exactamente la SC de tratamientos del DCA.
    """
    res = {}
    for org in df.microorganismo.unique():
        sub = df[df.microorganismo == org]
        m_dca = ols("recuento ~ C(trat)", data=sub).fit()
        t_dca = anova_lm(m_dca, typ=2)
        m_fac = ols("recuento ~ C(tiempo)*C(altura)", data=sub).fit()
        t_fac = anova_lm(m_fac, typ=2)
        res[org] = {
            "sc_trat_dca": float(t_dca.loc["C(trat)", "sum_sq"]),
            "sse_dca": float(t_dca.loc["Residual", "sum_sq"]),
            "r2_dca": float(m_dca.rsquared),
            "sc_tiempo": float(t_fac.loc["C(tiempo)", "sum_sq"]),
            "sc_altura": float(t_fac.loc["C(altura)", "sum_sq"]),
            "sc_inter": float(t_fac.loc["C(tiempo):C(altura)", "sum_sq"]),
            "sse_fac": float(t_fac.loc["Residual", "sum_sq"]),
            "r2_fac": float(m_fac.rsquared),
        }
    return res


def supuestos(df: pd.DataFrame, modelo) -> dict:
    """Normalidad (Shapiro-Wilk) y homogeneidad de varianzas (Levene)."""
    W, p_sh = stats.shapiro(modelo.resid)
    grupos = [g.recuento.values for _, g in df.groupby(["microorganismo", "trat"])]
    L, p_lev = stats.levene(*grupos, center="median")
    return {"shapiro_W": float(W), "shapiro_p": float(p_sh),
            "levene_stat": float(L), "levene_p": float(p_lev)}


def mejor_tratamiento(df: pd.DataFrame) -> dict:
    """Tratamiento de menor recuento medio por microorganismo."""
    out = {}
    for org in df.microorganismo.unique():
        medias = df[df.microorganismo == org].groupby("trat").recuento.mean().sort_values()
        out[org] = (medias.index[0], float(medias.iloc[0]))
    return out


def analisis_completo(path: Path | str = DATOS) -> dict:
    """Conveniencia: ejecuta todo el flujo y devuelve un diccionario de resultados."""
    df = cargar_datos(path)
    modelo, tabla = ajustar_factorial(df)
    SC, GL, MS = cuadrados_medios(tabla)
    return {
        "df": df, "modelo": modelo, "tabla": tabla,
        "SC": SC, "GL": GL, "MS": MS,
        "fijos": modelo_fijos(MS, GL),
        "aleatorios": modelo_aleatorios(MS, GL),
        "mixtos": modelo_mixtos(MS, GL),
        "componentes": componentes_varianza(MS),
        "dca": comparacion_dca(df),
        "supuestos": supuestos(df, modelo),
        "mejor": mejor_tratamiento(df),
        "r2": float(modelo.rsquared), "r2_adj": float(modelo.rsquared_adj),
    }


if __name__ == "__main__":  # resumen rápido por consola
    r = analisis_completo()
    print(f"N = {len(r['df'])} | R2 = {r['r2']:.4f}")
    print("Fijos  B×C:", round(r["fijos"]["BC"][0], 3))
    print("Mixtos B×C:", round(r["mixtos"]["BC"][0], 3))
    print("Comp. varianza BC:", round(r["componentes"]["BC"], 3))
    print("Mejor tratamiento:", r["mejor"])
