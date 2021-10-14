---
layout: default
title: Home
nav_order: 1
description: "Just the Docs is a responsive Jekyll theme with built-in search that is easily customizable and hosted on GitHub Pages."
permalink: /
---

# Pesquisa de equilíbrio de mercado em ambiente competitivo

Equipe de Modelos - PSR

## Introdução

Esta trabalho faz parte de uma pesquisa de mercado solicitada pelo ministério da economia para avaliar os melhores desenhos de mercado para o setor elétrico. Busca-se responder sobre a viabilidade de implantação de um mercado livre de intervenções governamentais no mercado de energia elétrica brasileiro. O modelo de otimização desenvolvido tenta descrever a dinâmica de expansão dos agentes geradores de energia em um mercado livre de intervenções. As interações entre os agentes e os consumidores se dão no âmbito do mercado de curto prazo.

## Mecanismo de competição

O mecanismo de competição construído pelo modelo parte da interação entre a garantia física total do sistema e o preço spot da energia elétrica. O objetivo ́e minimizar os custos de expansão da  geração, o  que  significa  construir  o  mínimo  necessário  para  poder  atender  a  demanda  durante um determinado período de tempo. Entretanto, o preço spot ́e decrescente com a garantia física. Uma menor garantia física total disponível implica num preço spot maior.   Um preço spot maior, por sua vez, atrai os investimentos para o setor da geração, induzindo crescimento da garantia física total. Um crescimento da garantia física total reduz o preço spot.  Percebe-se que há uma relação de compromisso entre a garantia física total e o preço da energia elétrica no mercado.  O modelo busca encontrar um ponto de equilíbrio neste sistema de forma a minimizar o preço da energia para o consumidor.

## Modelo desenvolvido

O modelo de otimização implementado no código é o seguinte:

<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

<div>
\[
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
\]
</div>

<div>
\[
\begin{aligned}
VPL_{i}(GF, Q_{i}) = \lambda_{g} E \left( \widehat{VPL}_{i}(GF, Q_{i}, c_{i}, g_{i, t, s}, \pi_{t, s})) \right) +
 (1 - \lambda_{g}) CVaR\left( \widehat{VPL}_{i}(GF, Q_{i}, c_{i}, g_{i, t, s}, \pi_{t, s}) \right)
 \end{aligned}
 \]
</div>

<div>
\[
\begin{aligned}
\widehat{VPL}_{i}(GF, Q_{i}, c_{i}, g_{i, t, s}, \pi_{t, s})) = \sum_{t} \left[g_{i, t, :} (\pi_{t, :} - c_{i}) + \left(\frac{1}{T}\sum_{t} \left[ E(\pi_{t, :}) \right] - \pi_{t, :}\right) Q_{i} \right]
\end{aligned}
\]
</div>

<div>
\[
\begin{aligned}
\Phi(GF, Q, d, \pi_{t, s}, \lambda_{d}, \alpha) = \lambda_{d} E \left( \widehat{\Phi}(GF, Q, d, \pi_{t, s})) \right) +
 (1 - \lambda_{d}) CVaR_\alpha\left(\widehat{\Phi}(GF, Q, d, \pi_{t, s})) \right)
 \end{aligned}
 \]
</div>

<div>
\[
\begin{aligned}
\widehat{\Phi}(GF, Q, d, \pi_{t, s}) = \sum_{t} \left[ (d-Q)\pi_{t, :} +  \frac{1}{T}\sum_{t} \left[ E(\pi_{t, :}) \right] Q\right]
\end{aligned}
\]
</div>

Onde:
* <div><nobr>\(d\): demanda de energia \(MW]\)</nobr></div>
* <div>\(GF\)</div>: garantia física total do sistema (somente usinas viabilizadas pela demanda) $[MWmédios]$
* <div>\(GF_i\)</div>: garantia física da usina $i$ $[MWmédios]$
* <div>\(Q\)</div>: montante total contratado pela demanda $[Wh]$
* <div>\(Q_i\)</div>: montante contratado pela demanda referente à usina $i$ $[Wh]$
* <div>\(x_i\)</div>: variável binária que indica se a demanda está disposta a viabilizar financeiramente à usina $i$
* <div>\(phi_i\)</div>: pagamento por fora feito à usina $i$
* M: constante aprox. $\infty$
* <div>\(VPL_i)</div>: valor presente líquido referente à usina $i$
* $\Phi$: gasto esperado pela demanda com a contratação das usinas 