# Problema
* Caixa bidimensional com N átomos
* Reflexão elástica nas paredes
* Distribuição espacial inicial aleatória
* Distribuição maxwelliana para as velocidades
    $$
    \rho(v) \propto  \exp\bigg(\frac{mv^2}{2k_bT}\bigg),~~~ m=1,~k_b=1,~v^2 = v_x^2 + v_y^2
    $$
* Potencial de interação entre as partículas
    $$
    \begin{aligned}
    V(r) &= \frac{1}{r^{2n}} - \frac{b}{r^n} \\
    f(r) &= -\frac{dV}{dr} = \frac{2n}{r^{2n+1}} - \frac{bn}{r^{n-1}}
    \end{aligned}
    $$

    A força é nula quando
    $$
    f(r^*) = 0 \Rightarrow r^* = \bigg(\frac{2}{b}\bigg)^{1/n}
    $$