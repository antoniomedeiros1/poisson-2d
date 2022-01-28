# Solucionador da Equação de Poisson 2D por Diferenças Finitas 

Programa implementado em *C++* para resolver a equação de Poisson 2D no domínio [0, 1] x [0, 1].

Foi utilizado o método das Diferenças Finitas para obter aproximações para as derivadas parciais
de segunda ordem, dividindo o domínio em N x M pontos discretos.

A discretização através desse método nos leva a um sistema linear, que pode ser resolvido considerando-se as 
condições de contorno u(0, y) = u(1, y) = u(x, 0) = 0, u(x, 1) = 1. O sistema linear foi resolvido utilizando métodos 
iterativos como Jacobi e SOR, considerando uma tolerância.

## Compilação e execução
Para compilar basta utilizar o comando:

``` make ```

Caso queira excluir os arquivos gerados pelo programa:

``` make clean ```

A execução do programa exige que seja passado como argumento a tolerância da solução:
``` ./main tol ```

Alguns shell scripts foram inclusos no repositório com objetivo te testar vários casos (variando a tolerancia e quantidade de threads) automaticamente.

## Resultados
Os resultados do programa, bem como uma análise do desempenho de cada um dos métodos pode ser encontrado no notebook *report.ipynb*.
