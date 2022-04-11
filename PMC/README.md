# Pesquisa de equilíbrio de mercado em ambiente competitivo

Equipe de Modelos - PSR

## Introdução

Esta trabalho faz parte de uma pesquisa de mercado solicitada pelo ministério da economia para avaliar os melhores desenhos de mercado para o setor elétrico. Busca-se responder sobre a viabilidade de implantação de um mercado livre de intervenções governamentais no mercado de energia elétrica brasileiro. O modelo de otimização desenvolvido tenta descrever a dinâmica de comercialização da energia em um mercado livre de intervenções. As interações entre os agentes e os consumidores se dão no âmbito do mercado de curto prazo.

## Modelo desenvolvido

O modelo de otimização implementado no código é o seguinte:

$$
\begin{aligned}
    & \underset{GF, Q, Q_i, x, \phi}{\text{min:}}
    & & \Phi \left(GF, Q \right) + \sum_{i=1}^{N} \phi_i\\
    & \text{sujeita à} & & GF - \sum_{i=1}^{N} GF_i x_i = 0\\
    & & &  VPL_i(GF, Q_i) + \phi_i \geq -M\left(1-x_i\right) \; \;(i=1, 2, ..., N)\\
    & & & Q - \sum_{i=1}^{N} Q_i= 0\\
    & & & Q_i \leq d x_i \\
    & & & \phi_i \geq 0 \; \;(i=1, 2, ..., N)\\
    & & & Q_i \geq 0 \; \;(i=1, 2, ..., N)\\
\end{aligned}
$$

$$
VPL_{i}(GF, Q_{i}) = \lambda_{g} E \left( \widehat{VPL}_{i}(GF, Q_{i}, c_{i}, g_{i, t, s}, \pi_{t, s})) \right) +
 (1 - \lambda_{g}) CVaR\left( \widehat{VPL}_{i}(GF, Q_{i}, c_{i}, g_{i, t, s}, \pi_{t, s}) \right)
\\

\widehat{VPL}_{i}(GF, Q_{i}, c_{i}, g_{i, t, s}, \pi_{t, s})) = \sum_{t} \left[g_{i, t, :} (\pi_{t, :} - c_{i}) + \left(\frac{1}{T}\sum_{t} \left[ E(\pi_{t, :}) \right] - \pi_{t, :}\right) Q_{i} \right]
\\

\Phi(GF, Q, d, \pi_{t, s}, \lambda_{d}, \alpha) = \lambda_{d} E \left( \widehat{\Phi}(GF, Q, d, \pi_{t, s})) \right) +
 (1 - \lambda_{d}) CVaR_\alpha\left(\widehat{\Phi}(GF, Q, d, \pi_{t, s})) \right)
\\

\widehat{\Phi}(GF, Q, d, \pi_{t, s}) = \sum_{t} \left[ (d-Q)\pi_{t, :} +  \frac{1}{T}\sum_{t} \left[ E(\pi_{t, :}) \right] Q\right]
$$

Onde:
* $d$: demanda de energia $[MW]$
* $GF$: garantia física total do sistema (somente usinas viabilizadas pela demanda) $[MWmédios]$
* $GF_i$: garantia física da usina $i$ $[MWmédios]$
* $Q$: montante total contratado pela demanda $[Wh]$
* $Q_i$: montante contratado pela demanda referente à usina $i$ $[Wh]$
* $x_i$: variável binária que indica se a demanda está disposta a viabilizar financeiramente à usina $i$
* $\phi_i$: pagamento por fora feito à usina $i$
* M: constante aprox. $\infty$
* $VPL_i$: valor presente líquido referente à usina $i$
* $\Phi$: gasto esperado pela demanda com a contratação das usinas