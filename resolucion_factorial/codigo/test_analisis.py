# -*- coding: utf-8 -*-
"""Pruebas unitarias del núcleo de cálculo (analisis.py).

Ejecutar desde la carpeta ``codigo``:
    pytest -v

Cada prueba fija un valor verificado contra los datos crudos del primer examen,
de modo que cualquier cambio accidental en una fórmula (p. ej. un denominador de
EMS o un componente de varianza) se detecte de inmediato.
"""
import math

import pytest

import analisis as an

TOL = 1e-3  # tolerancia para comparaciones numéricas


# --------------------------------------------------------------------------- #
# Fixtures: se ajusta el modelo una sola vez y se comparte entre pruebas
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def datos():
    return an.cargar_datos()


@pytest.fixture(scope="module")
def ajuste(datos):
    modelo, tabla = an.ajustar_factorial(datos)
    SC, GL, MS = an.cuadrados_medios(tabla)
    return {"modelo": modelo, "tabla": tabla, "SC": SC, "GL": GL, "MS": MS}


# --------------------------------------------------------------------------- #
# 1. Estructura del diseño
# --------------------------------------------------------------------------- #
def test_numero_de_observaciones(datos):
    assert len(datos) == an.A_NIV * an.B_NIV * an.C_NIV * an.N_REP == 72


def test_diseno_balanceado(datos):
    conteos = datos.groupby(["microorganismo", "tiempo", "altura"]).size()
    assert conteos.nunique() == 1
    assert conteos.iloc[0] == an.N_REP == 4


def test_niveles_de_factores(datos):
    assert datos.microorganismo.nunique() == 2
    assert sorted(datos.tiempo.unique()) == [1, 3, 5]
    assert sorted(datos.altura.unique()) == [5, 10, 20]


# --------------------------------------------------------------------------- #
# 2. Grados de libertad y cuadrados medios
# --------------------------------------------------------------------------- #
def test_grados_de_libertad(ajuste):
    esperado = {"A": 1, "B": 2, "C": 2, "AB": 2, "AC": 2, "BC": 4, "ABC": 4, "E": 54}
    assert ajuste["GL"] == esperado


def test_cuadrados_medios(ajuste):
    MS = ajuste["MS"]
    esperado = {"A": 19.0139, "B": 86.0972, "C": 58.6806, "AB": 22.2639,
                "AC": 0.8472, "BC": 94.5347, "ABC": 2.4097, "E": 1.7917}
    for k, v in esperado.items():
        assert MS[k] == pytest.approx(v, abs=1e-3)


def test_suma_de_cuadrados_total(ajuste):
    SC = ajuste["SC"]
    total = sum(SC.values())
    assert total == pytest.approx(839.3194, abs=1e-2)


# --------------------------------------------------------------------------- #
# 3. Modelo I — efectos fijos (todo vs. error)
# --------------------------------------------------------------------------- #
def test_modelo_fijos_valores_F(ajuste):
    f = an.modelo_fijos(ajuste["MS"], ajuste["GL"])
    esperado_F = {"A": 10.612, "B": 48.054, "C": 32.752, "AB": 12.426,
                  "AC": 0.473, "BC": 52.764, "ABC": 1.345}
    for k, F in esperado_F.items():
        assert f[k][0] == pytest.approx(esperado_F[k], abs=1e-2)
        assert f[k][2] == 54  # denominador = gl del error


def test_modelo_fijos_significancia(ajuste):
    f = an.modelo_fijos(ajuste["MS"], ajuste["GL"])
    # Significativos: A, B, C, AB, BC ; no significativos: AC, ABC
    assert f["BC"][3] < 0.001
    assert all(f[k][3] < 0.05 for k in ["A", "B", "C", "AB", "BC"])
    assert all(f[k][3] > 0.05 for k in ["AC", "ABC"])


# --------------------------------------------------------------------------- #
# 4. Modelo II — efectos aleatorios (denominadores EMS + Satterthwaite)
# --------------------------------------------------------------------------- #
def test_modelo_aleatorios_denominadores(ajuste):
    MS, GL = ajuste["MS"], ajuste["GL"]
    al = an.modelo_aleatorios(MS, GL)
    # Interacciones dobles contra la triple
    assert al["AB"][0] == pytest.approx(MS["AB"] / MS["ABC"], abs=TOL)
    assert al["AB"][0] == pytest.approx(9.239, abs=1e-2)
    assert al["BC"][0] == pytest.approx(39.231, abs=1e-2)
    assert al["BC"][2] == 4  # gl denominador = gl de ABC
    # Triple contra el error
    assert al["ABC"][0] == pytest.approx(1.345, abs=1e-2)
    assert al["ABC"][2] == 54


def test_modelo_aleatorios_satterthwaite(ajuste):
    al = an.modelo_aleatorios(ajuste["MS"], ajuste["GL"])
    # Efectos principales: F pequeño y NO significativo bajo el modelo aleatorio
    assert al["A"][0] == pytest.approx(0.918, abs=1e-2)
    assert al["A"][2] == pytest.approx(1.72, abs=0.05)  # gl Satterthwaite
    for k in ["A", "B", "C"]:
        assert al[k][3] > 0.05


# --------------------------------------------------------------------------- #
# 5. Modelo III — efectos mixtos (A fijo; B,C aleatorios; restringido)
# --------------------------------------------------------------------------- #
def test_modelo_mixtos_denominadores(ajuste):
    MS, GL = ajuste["MS"], ajuste["GL"]
    mx = an.modelo_mixtos(MS, GL)
    # B y C contra CM_BC
    assert mx["B"][0] == pytest.approx(MS["B"] / MS["BC"], abs=TOL)
    assert mx["B"][0] == pytest.approx(0.911, abs=1e-2)
    assert mx["C"][0] == pytest.approx(0.621, abs=1e-2)
    assert mx["B"][2] == 4 and mx["C"][2] == 4
    # BC contra el error -> sigue siendo muy significativa
    assert mx["BC"][0] == pytest.approx(52.764, abs=1e-2)
    assert mx["BC"][2] == 54
    assert mx["BC"][3] < 0.001


def test_interaccion_BC_robusta_en_los_tres_modelos(ajuste):
    """La interacción tiempo×altura es significativa en fijos, aleatorios y mixtos."""
    MS, GL = ajuste["MS"], ajuste["GL"]
    f = an.modelo_fijos(MS, GL)
    al = an.modelo_aleatorios(MS, GL)
    mx = an.modelo_mixtos(MS, GL)
    assert f["BC"][3] < 0.05
    assert al["BC"][3] < 0.05
    assert mx["BC"][3] < 0.05


def test_efectos_principales_pierden_significancia(ajuste):
    """Efectos principales: significativos en fijos, no en aleatorios/mixtos."""
    MS, GL = ajuste["MS"], ajuste["GL"]
    f = an.modelo_fijos(MS, GL)
    al = an.modelo_aleatorios(MS, GL)
    mx = an.modelo_mixtos(MS, GL)
    for k in ["B", "C"]:
        assert f[k][3] < 0.05
        assert al[k][3] > 0.05
        assert mx[k][3] > 0.05


# --------------------------------------------------------------------------- #
# 6. Componentes de varianza (modelo aleatorio)
# --------------------------------------------------------------------------- #
def test_componentes_varianza(ajuste):
    vc = an.componentes_varianza(ajuste["MS"])
    assert vc["BC"] == pytest.approx(11.5156, abs=1e-3)   # domina
    assert vc["AB"] == pytest.approx(1.6545, abs=1e-3)
    assert vc["error"] == pytest.approx(1.7917, abs=1e-3)
    assert vc["ABC"] == pytest.approx(0.1545, abs=1e-3)


def test_componentes_negativos(ajuste):
    """Varios componentes son negativos => el supuesto aleatorio no es coherente."""
    vc = an.componentes_varianza(ajuste["MS"])
    for k in ["A", "B", "C", "AC"]:
        assert vc[k] < 0


def test_BC_domina_la_varianza(ajuste):
    vc = an.componentes_varianza(ajuste["MS"])
    pos = {k: max(v, 0) for k, v in vc.items()}
    frac_bc = pos["BC"] / sum(pos.values())
    assert frac_bc > 0.70  # ~76 %


# --------------------------------------------------------------------------- #
# 7. Comparación factorial vs. DCA del primer examen
# --------------------------------------------------------------------------- #
def test_factorial_descompone_la_SC_del_DCA(datos):
    dca = an.comparacion_dca(datos)
    for org, r in dca.items():
        suma = r["sc_tiempo"] + r["sc_altura"] + r["sc_inter"]
        assert suma == pytest.approx(r["sc_trat_dca"], abs=1e-6)


def test_factorial_conserva_SSE_y_R2(datos):
    dca = an.comparacion_dca(datos)
    for org, r in dca.items():
        assert r["sse_fac"] == pytest.approx(r["sse_dca"], abs=1e-6)
        assert r["r2_fac"] == pytest.approx(r["r2_dca"], abs=1e-9)


def test_valores_DCA_conocidos(datos):
    dca = an.comparacion_dca(datos)
    assert dca["E. coli"]["sc_trat_dca"] == pytest.approx(267.50, abs=1e-2)
    assert dca["S. typhimurium"]["sc_trat_dca"] == pytest.approx(456.06, abs=1e-2)


# --------------------------------------------------------------------------- #
# 8. Supuestos y tratamiento óptimo
# --------------------------------------------------------------------------- #
def test_supuestos(ajuste, datos):
    s = an.supuestos(datos, ajuste["modelo"])
    assert s["shapiro_p"] == pytest.approx(0.0685, abs=2e-3)
    assert s["levene_p"] == pytest.approx(0.9146, abs=2e-3)
    assert s["shapiro_p"] > 0.05   # normalidad global no se rechaza
    assert s["levene_p"] > 0.05    # homogeneidad de varianzas


def test_mejor_tratamiento(datos):
    mejor = an.mejor_tratamiento(datos)
    assert mejor["E. coli"][0] == "T3_H10"
    assert mejor["E. coli"][1] == pytest.approx(6.75, abs=1e-9)
    assert mejor["S. typhimurium"][0] == "T3_H10"
    assert mejor["S. typhimurium"][1] == pytest.approx(6.00, abs=1e-9)


# --------------------------------------------------------------------------- #
# 9. Consistencia interna de las pruebas F (mismo numerador entre modelos)
# --------------------------------------------------------------------------- #
def test_consistencia_F_ABC_entre_modelos(ajuste):
    """ABC se prueba siempre contra el error: mismo F en los tres modelos."""
    MS, GL = ajuste["MS"], ajuste["GL"]
    f = an.modelo_fijos(MS, GL)
    al = an.modelo_aleatorios(MS, GL)
    mx = an.modelo_mixtos(MS, GL)
    assert f["ABC"][0] == pytest.approx(al["ABC"][0], abs=TOL)
    assert f["ABC"][0] == pytest.approx(mx["ABC"][0], abs=TOL)
