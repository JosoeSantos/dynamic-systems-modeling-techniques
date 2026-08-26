import math

import matplotlib.pyplot as plt
import numpy as np

from functools import partial

# plt.rcParams["text.usetex"] = True
# ---------------------------------------------------------------------------
# 1. FIXED (process/coffee-related) PARAMETERS  -- Table 5, left column
# ---------------------------------------------------------------------------
FIXED_PARAMS = {
    "A": 116200.0,  # Arrhenius pre-factor for exothermic reaction   [kJ/kg]
    "cm": 0.418,  # specific heat of the metal                    [kJ/(kg K)]
    "cw": 5.0,  # partial specific heat of water in the bean    [kJ/(kg K)]
    "Db": 7.65e-3,  # bean diameter                                 [m]
    "Ha_R": 5500.0,  # activation energy / gas constant (Ha/R)       [K]
    "Het": 232.0,  # total reaction heat                           [kJ/kg]
    "k1": 4.32e-9,  # Arrhenius-type moisture-loss parameter
    "k2": 9889.0,  # Arrhenius-type moisture-loss parameter        [K]
    "Kt": 0.01,  # bean-temperature sensor time constant         [1/s]
    "mb": 1.5e-4,  # mass of a single bean                         [kg]
    "lam": 2790.0,  # latent heat of vaporization of bean moisture  [kJ/kg]
}

# ---------------------------------------------------------------------------
# 2. SCALABLE (geometry-related) PARAMETERS  -- Table 5, right column
# ---------------------------------------------------------------------------
PLANTS = {
    "120kg": {
        "Mb": 120.0,
        "Dch": 1.24,
        "Hflap": 0.3,
        "Lch": 1.335,
        "Mm": 2000.0,
        "Sflap": 0.1,
    },
    "360kg": {
        "Mb": 360.0,
        "Dch": 1.90,
        "Hflap": 0.3,
        "Lch": 2.04,
        "Mm": 7000.0,
        "Sflap": 0.1,
    },
}


# ---------------------------------------------------------------------------
# 3. IDENTIFIED PARAMETERS -- Table 1 (identified once on the 120 kg plant,
#    then reused unchanged, i.e. "transplanted", on the 360 kg plant)
# ---------------------------------------------------------------------------
IDENTIFIED_PARAMS = {
    "hgm": 0.0100,  # gas -> metal heat transfer coefficient   [kW/(m^2 K)] (see UNIT NOTE)
    "hbm": 0.0254,  # metal -> bean heat transfer coefficient  [kW/(m^2 K)]
    "Pbm": 0.5793,  # fraction of a bean's surface in contact with the metal [-]
}

# ---------------------------------------------------------------------------
# 4. Air specific-heat polynomial, Eq. (2)  (originally J/(kg K), Tgi in degC)
# ---------------------------------------------------------------------------
_ALPHA = [
    1.0839e3,
    -7.2075e-1,
    2.1034e-3,
    -2.3267e-6,
    1.3621e-9,
    -4.1550e-13,
    5.3091e-17,
]


def he_coeff(X):
    """Gas -> bean heat transfer coefficient, Eq. (3), function of moisture X."""
    return 0.49 - 0.443 * np.exp(-0.206 * X)


def cg_air(Tgi_degC):
    """Specific heat capacity of drying air Eq. (2). kJ/(kg K)"""
    Tk = Tgi_degC + 273.15
    cg_J = sum(a * Tk**i for i, a in enumerate(_ALPHA))
    return cg_J / 1000.0


def geometric_params(plant):
    """Geometric params based in the plant format and beans geometric params"""
    Pbm, Mb = IDENTIFIED_PARAMS["Pbm"], plant["Mb"]
    mb, Db = FIXED_PARAMS["mb"], FIXED_PARAMS["Db"]
    # Area of contact gas-metal
    Agm = (
        math.pi
        * plant["Dch"]
        * (
            plant["Lch"]
            + plant["Hflap"] * plant["Lch"] / plant["Sflap"]
            + plant["Dch"] / 2
        )
    )
    # Total beans surface area
    Ab = (Mb / mb) * math.pi * Db**2
    # Gas to bean heat transfer Area
    Agb = Ab * (1 - Pbm)
    # Area of contact bean metal
    Abm = Ab * Pbm

    return (Agm, Ab, Agb, Abm)


def compute_F(hgm, he, Agm, Agb):
    """Ration of gas-metal and gas-beans thermal resistances"""
    return hgm * Agm / (he * Agb)


def temp_gas_outlet(Tgi, Tb, F, Tm, he, Agb, Gg, cg):
    Gg = max(Gg, 1e-9)  # guard against div-by-zero, not just negative
    driving = (Tgi - (Tb + F * Tm) / (1 + F)) 
    coeff = math.exp(-he * Agb  / (Gg * cg))
    return (Tgi - driving) * (1 - coeff)


def dynamic_model_derivatives(x, u, t, fixed, identified, geometric):
    (
        Tb,
        Tm,
        X,
        He,
        Ta,
    ) = x
    Gg, Tgi = u
    Agm, Ab, Agb, Abm = geometric_params(geometric)
    hgm = identified["hgm"]
    hbm = identified["hbm"]
    Pbm = identified["Pbm"]

    k1 = fixed["k1"]
    k2 = fixed["k2"]
    Db = fixed["Db"]
    A = fixed["A"]
    Het = fixed["Het"]
    Ha_R = fixed["Ha_R"]
    Kt = fixed["Kt"]

    he = he_coeff(X)
    cg = cg_air(Tgi)
    F = compute_F(hgm, he, Agm, Agb)

    Tgo = temp_gas_outlet(Tgi, Tb, F, Tm, he, Agb, Gg, cg)

    # Dynamic updates

    # Gas bean heat transfer Eq. 9
    Qgb = Gg * cg * (Tgi - Tgo)

    # Gas metal heat rate Eq. 10
    Qgm = F * (he * Agb * (Tb - Tm) + Qgb) / (1 + F)

    # Bean metal heat transfer rate Eq. 11
    Qbm = hbm * Abm * (Tm - Tb)

    ## Mass of beans
    Mb = geometric["Mb"]

    X0 = geometric["X0"]

    ## Dry beans mass


    Qe = 0  # Negligible heat loss

    cs = 1.099 + 0.007 * Tb
    cb = (cs + fixed["cw"] * X) / (1 + X)

    Qr = A * ((Het - He) / Het) * math.exp(-Ha_R / (Tb + 273.15))

    # x_3_dot
    # Moisture
    X_dot = -(k1 / Db**2) * math.exp(-k2 / (Tb + 273.15))
    Mbd = Mb / (1 + X)
    # x_1_dot
    Tb_dot = (Qgb - Qgm + Qbm + Mbd * (Qr + fixed["lam"] * X_dot)) / (
        Mbd * (1 + X) * cb
    )

    # x_2_dot
    # Metal temp Eq. 12
    Tm_dot = (Qgm - Qbm + Qe) / (geometric["Mm"] * fixed["cm"])

    # x_4_dot
    # Amount of heat produced by kilogram of coffee thus far
    He_dot = Qr

    Ta_dot = Kt * (Tb - Ta)

    return np.array([Tb_dot, Tm_dot, X_dot, He_dot, Ta_dot])


## Implementar sim
## IMPLEMENTAR RK4


def rk4(dydt, x0, u, t, h):
    """Function Runge Kunta of 4th order

    dydt derivative of system
    x0 last point
    u system inputs
    h integration step single
    t starting time

    """
    k1 = dydt(x0, u, t)

    k2 = dydt(x0 + 0.5 * h * k1, u, t + 0.5 * h)

    k3 = dydt(x0 + 0.5 * h * k2, u, t + 0.5 * h)

    k4 = dydt(x0 + h * k3, u, t + h)

    return x0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def constant_u(t, temp=120):
    return np.array([np.full(len(t), 0.2), np.full(len(t), temp)])


def step_u(t):
    Tgi = np.full(len(t), 120)
    Tgi[math.floor(len(t) / 2) : len(t)] = 180
    return np.array([np.full(len(t), 0.2), Tgi])


def impulse_u(t, dt_len):
    pulse_width = 1
    pulse_time = math.floor(len(t) / 2)

    Tgi = np.full(len(t), 120)
    Tgi[pulse_time : pulse_time + pulse_width] = 1 / dt_len

    return np.array([np.full(len(t), 0.2), Tgi])


def wave_u(t):
    # Create 100 points from 0 to 2*pi (one full sine wave cycle)
    x = np.linspace(0, 2 * np.pi, len(t))
    Tgi = np.sin(x) * 30 + 90
    return np.array([np.full(len(t), 0.2), Tgi])


def random_u(t):
    rng = np.random.default_rng()
    Tgi = rng.random(len(t)) * 30 + 90
    return np.array([np.full(len(t), 0.2), Tgi])


# ---------------------------------------------------------------------------
# NÚCLEO DE SIMULAÇÃO (separado da plotagem, para permitir reaproveitamento
# em comparações, como no teste do princípio da superposição)
# ---------------------------------------------------------------------------
def run_simulation(time, plant, uf, h, x0):
    """Executa a simulação e retorna (t, x, y), sem plotar nada.

    time = { t_start, t_end }
    plant = "120kg" | "360kg"
    uf = função geradora de entrada u(t)
    h = passo de integração
    x0 = condição inicial
    """
    t = np.arange(time["t_start"], time["t_end"], h)

    x = np.zeros((len(x0), len(t)))
    y = np.zeros(len(t))
    u = uf(t)

    # cópia para não mutar o dicionário global PLANTS entre chamadas
    geo = dict(PLANTS[plant])
    geo["X0"] = x0[2]

    dydt = partial(
        dynamic_model_derivatives,
        fixed=FIXED_PARAMS,
        geometric=geo,
        identified=IDENTIFIED_PARAMS,
    )

    x[:, 0] = x0
    for i in range(1, len(t)):
        x[:, i] = rk4(dydt, x[:, i - 1], u[:, i], t[i - 1], h)

        (Tb, Tm, X, He, Ta) = x[:, i]
        Gg, Tgi = u[:, i]
        hgm = IDENTIFIED_PARAMS["hgm"]
        Agm, Ab, Agb, Abm = geometric_params(geo)
        he = he_coeff(X)
        cg = cg_air(Tgi)
        F = compute_F(hgm, he, Agm, Agb)
        y[i] = temp_gas_outlet(Tgi, Tb, F, Tm, he, Agb, Gg, cg)

    return t, x, y, u


def simulate(
    time,
    plant,
    uf,
    h,
    x0,
    title="Câmara de torrefação de café",
):
    """Executa a simulação e plota os resultados (comportamento original)."""

    t, x, y, u = run_simulation(time, plant, uf, h, x0)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle(title)
    axes[0, 0].plot(t, y, label="Tgo")
    axes[0, 0].plot(t, u[1, :], label="Desired Tgi")
    axes[0, 0].set_title("Temperature GO (output)")
    axes[0, 0].set_xlabel("Time [s]")
    axes[0, 0].set_ylabel("Temperature [%]")
    axes[0, 0].legend()

    axes[0, 1].plot(t, x[1, :], label="Tm")
    axes[0, 1].plot(t, x[4, :], label="Ta")
    axes[0, 1].set_title("Temperature machine (estimated)")
    axes[0, 1].set_xlabel("Time [s]")
    axes[0, 1].set_ylabel("Temperature [%]")
    axes[0, 1].legend()

    axes[1, 0].plot(t, x[3, :], label="He")
    axes[1, 0].set_title("Heat produced per Kg of coffee")
    axes[1, 0].set_xlabel("Time [s]")
    axes[1, 0].set_ylabel("Temperature [%]")
    axes[1, 0].legend()

    axes[1, 1].plot(t, x[2, :], label="Moisture")
    axes[1, 1].set_title("Moisture loss over time")
    axes[1, 1].set_xlabel("Time [s]")
    axes[1, 1].set_ylabel("Moisture content [%]")
    axes[1, 1].legend()
    fig.tight_layout()
    plt.show()

    return t, x, y, u


# ---------------------------------------------------------------------------
# TESTE DO PRINCÍPIO DA SUPERPOSIÇÃO
#
# Um sistema é LINEAR se, e somente se, satisfizer simultaneamente:
#   (i)  Aditividade:      y(u1 + u2) = y(u1) + y(u2)
#   (ii) Homogeneidade:    y(k*u)     = k*y(u)
#
# Como o modelo tem um ponto de operação (equilíbrio) diferente de zero,
# as duas propriedades são testadas em variáveis de DESVIO em torno de uma
# entrada de referência (baseline) constante, partindo sempre da MESMA
# condição inicial x0:
#
#   (i)  [y(u0+d1) - y(u0)] + [y(u0+d2) - y(u0)]  vs.  y(u0+d1+d2) - y(u0)
#   (ii) k * [y(u0+d) - y(u0)]                    vs.  y(u0+k*d) - y(u0)
#
# Se o sistema fosse linear, as duas curvas de cada teste coincidiriam
# perfeitamente (erro = 0). Qualquer diferença comprova a NÃO LINEARIDADE.
# ---------------------------------------------------------------------------
def _const_input_factory(Gg_const):
    def u_const(t, Tgi):
        return np.array([np.full(len(t), Gg_const), np.full(len(t), Tgi)])

    return u_const


def test_additivity(
    time,
    plant,
    h,
    x0,
    Gg_const=0.2,
    Tgi_base=120,
    delta1=30,
    delta2=15,
    title="Teste de Aditividade (Superposição)",
):
    """Testa y(u1+u2) = y(u1) + y(u2) (em variáveis de desvio)."""

    u_const = _const_input_factory(Gg_const)

    u0 = partial(u_const, Tgi=Tgi_base)
    u1 = partial(u_const, Tgi=Tgi_base + delta1)
    u2 = partial(u_const, Tgi=Tgi_base + delta2)
    u12 = partial(u_const, Tgi=Tgi_base + delta1 + delta2)

    t, _, y0, _ = run_simulation(time, plant, u0, h, x0)
    _, _, y1, _ = run_simulation(time, plant, u1, h, x0)
    _, _, y2, _ = run_simulation(time, plant, u2, h, x0)
    _, _, y12, _ = run_simulation(time, plant, u12, h, x0)

    # resposta que a superposição PREVERIA, caso o sistema fosse linear
    y_previsto = y1 + y2 - y0
    erro = y12 - y_previsto

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(title)

    axes[0].plot(t, y12, label=r"$y(u_1+u_2)$  (resposta real do sistema)")
    axes[0].plot(
        t,
        y_previsto,
        "--",
        label=r"$y(u_1)+y(u_2)-y_0$  (previsto pela superposição)",
    )
    axes[0].set_title("Resposta real x resposta prevista pela superposição")
    axes[0].set_xlabel("Tempo [s]")
    axes[0].set_ylabel("Tgo [°C]")
    axes[0].legend()

    axes[1].plot(t, erro, color="tab:red")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_title("Erro de superposição (real - previsto)")
    axes[1].set_xlabel("Tempo [s]")
    axes[1].set_ylabel("Erro [°C]")

    fig.tight_layout()
    plt.show()

    erro_max = np.max(np.abs(erro))
    print(f"[Aditividade] Erro máximo de superposição: {erro_max:.4f} °C")

    return t, y0, y1, y2, y12, erro


def test_homogeneity(
    time,
    plant,
    h,
    x0,
    Gg_const=0.2,
    Tgi_base=120,
    delta=20,
    scale=2,
    title="Teste de Homogeneidade (Superposição)",
):
    """Testa y(k*u) = k*y(u) (em variáveis de desvio)."""

    u_const = _const_input_factory(Gg_const)

    u0 = partial(u_const, Tgi=Tgi_base)
    u1 = partial(u_const, Tgi=Tgi_base + delta)
    u_scaled = partial(u_const, Tgi=Tgi_base + scale * delta)

    t, _, y0, _ = run_simulation(time, plant, u0, h, x0)
    _, _, y1, _ = run_simulation(time, plant, u1, h, x0)
    _, _, y_scaled, _ = run_simulation(time, plant, u_scaled, h, x0)

    y_previsto = scale * (y1 - y0) + y0
    erro = y_scaled - y_previsto

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(title)

    axes[0].plot(t, y_scaled, label=rf"$y({scale}\cdot u)$  (resposta real)")
    axes[0].plot(
        t,
        y_previsto,
        "--",
        label=rf"${scale}\cdot[y(u)-y_0]+y_0$  (previsto)",
    )
    axes[0].set_title("Resposta real x resposta prevista pela homogeneidade")
    axes[0].set_xlabel("Tempo [s]")
    axes[0].set_ylabel("Tgo [°C]")
    axes[0].legend()

    axes[1].plot(t, erro, color="tab:red")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_title("Erro de homogeneidade (real - previsto)")
    axes[1].set_xlabel("Tempo [s]")
    axes[1].set_ylabel("Erro [°C]")

    fig.tight_layout()
    plt.show()

    erro_max = np.max(np.abs(erro))
    print(f"[Homogeneidade] Erro máximo de superposição: {erro_max:.4f} °C")

    return t, y0, y1, y_scaled, erro


# %% Op point and sim constants

X0 = 0.1

# Initial conditions from paper
# Tb, Tm, X, He, Ta,
x0 = np.array([30.0, 121.0, X0, 0.0, 30.0])
time_limits = {"t_start": 0, "t_end": 800}

integration_step = 1


# %%

impulse = partial(impulse_u, dt_len=integration_step / 800)


simulate(
    time_limits,
    "120kg",
    step_u,
    integration_step,
    x0,
    r"Câmara de torrefação de café - Função degrau",
)

simulate(
    time_limits,
    "120kg",
    impulse,
    integration_step,
    x0,
    r"Câmara de torrefação de café - Impulso",
)

simulate(
    time_limits,
    "120kg",
    wave_u,
    integration_step,
    x0,
    r"Câmara de torrefação de café - Entrada senoidal",
)

simulate(
    time_limits,
    "120kg",
    random_u,
    integration_step,
    x0,
    r"Câmara de torrefação de café - Entrada aleatória (ruído)",
)

# %% Princípio da superposição - visão inicial (entradas constantes isoladas)

simulate(
    time_limits,
    "120kg",
    constant_u,
    integration_step,
    x0,
    "Câmara de torrefação de café - Entrada constante 120",
)

constant_u_60 = partial(constant_u, temp=60)

simulate(
    time_limits,
    "120kg",
    constant_u_60,
    integration_step,
    x0,
    "Câmara de torrefação de café - Entrada constante 60",
)

# %% Princípio da superposição - teste formal (aditividade + homogeneidade)
#
# Aqui comparamos a resposta real do sistema não linear com a resposta que
# SERIA obtida se ele fosse linear (soma/escala das respostas individuais).
# A diferença entre as curvas (erro != 0) é a prova de que o sistema
# NÃO É LINEAR.

test_additivity(
    time_limits,
    "120kg",
    integration_step,
    x0,
    Gg_const=0.2,
    Tgi_base=120,
    delta1=60,
    delta2=60,
    title="Câmara de torrefação de café - Teste de Aditividade",
)

test_homogeneity(
    time_limits,
    "120kg",
    integration_step,
    x0,
    Gg_const=0.2,
    Tgi_base=120,
    delta=20,
    scale=2,
    title="Câmara de torrefação de café - Teste de Homogeneidade",
)
