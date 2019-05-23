# Computing with real numbers (in python)

The implementation is based on [Computing with Real Numbers][1]. Here's a quick breakdown for
the interested reader:

As our primitive, we model only numbers in the interval $[-1, 1]$. Other numbers will be represented
by scaling appropriately. We will work in base $2$, although this can be tweaked in the implementation.
Any power of two should work fine out of the box. A number for us will be an infinite stream of digits.

Consider an ordinary digit stream of base $2$, meaning you either get a $0$ or a $1$. The normal arithmetic
tells you that seeing a $0$ as the first digit means the rest of the number must be in the interval $[0, 0.5)$,
seeing a $1$ similarly restricts to $[0.5, 1)$. The problem with this approach is obvious once you try
to add two such numbers. When one of the numbers is a stream of $1$s and the other one is a stream of $0$ you
still can not possibly say if the result must begin with a $1$! Any $1$ in the second stream would lead to
a carry and the first digit of the result would be a $1$. If there is no $1$ in the second stream, there would be
no carry, and the first digit of the result would be a $0$.

As a remedy, we will allow *negative* digits, also, which we write as $\bar{1}$.
Formally, we will interpret a stream of digits as a so called [Linear fractional transformation][2], say $L(x)$,
and after seeing a finite number of digits, the real number of the whole stream is somewhere in the interval
$L([-1, 1])$. To be more precise, we will interpret a $\bar{1}$ as the transformation
$M_{\bar{1}}(x) = \frac{x - 1}{2}$, $0$ is interpreted as $M_0(x) = \frac{x}{2}$ and $1$ as $M_1(x) = \frac{x + 1}{2}$.
Seeing multiple digits after another $d e$ is simply the composition of those functions $M_{d e}(x) = (M_d \circ M_e)(x)$.

## Calculation of compositions of linear fractional transformations

Before we go to calculating anything, there is a nice trick for calculating with LFTs. Namely, if you have
two LFTs $A(x) = \frac{ax + b}{cx + d}$ and $B(x) = \frac{ux + v}{wx + y}$, then their composition can be calculated
easily. Suggestively reading each LFT as a matrix, their composition is just the normal matrix product!

$$
\begin{aligned}

A & \equiv \bigl( \begin{matrix} a & b \\ c & d \end{matrix} \bigr) \\

B & \equiv \bigl( \begin{matrix} u & v \\ w & y \end{matrix} \bigr) \\

(A \circ B) & \equiv \bigl( \begin{matrix} a & b \\ c & d \end{matrix} \bigr) * \bigl( \begin{matrix} u & v \\ w & y \end{matrix} \bigr)
            & = \bigl( \begin{matrix} a * u + b * w & a * v + b * y \\ c * u + d * w & c * v + d * y \end{matrix} \bigr)
\end{aligned}
$$

## A simple example of a single number

Let's do an example to get warmed up: Say, we observe the stream of digits $[0, \bar{1}, 1]$, then the number will
be contained somewhere in the following interval:

$$
\begin{aligned}
M_{0 \bar{1} 1} & = (M_0 \circ M_{\bar{1}} \circ M_1) \\
                & \equiv \bigl( \begin{matrix} 1 & 0 \\ 0 & 2 \end{matrix} \bigr)
                        * \bigl( \begin{matrix} 1 & -1 \\ 0 & 2 \end{matrix} \bigr)
                        * \bigl( \begin{matrix} 1 & 1 \\ 0 & 2 \end{matrix} \bigr) \\
                & = \bigl( \begin{matrix} 1 & -1 \\ 0 & 4 \end{matrix} \bigr)
                        * \bigl( \begin{matrix} 1 & 1 \\ 0 & 2 \end{matrix} \bigr) \\
                & = \bigl( \begin{matrix} 1 & -1 \\ 0 & 8 \end{matrix} \bigr) \\
                & \equiv (\frac{x - 1}{8}).
\end{aligned}
$$

We conclude that the number starting with $[0, \bar{1}, 1]$ must be somewhere in the interval $[-0.25, 0]$.

## Unary and binary operations

Since we now know how to interpret a single number, we can see if we can build up some operations on them. It would
be convenient, if those were also LFTs and indeed, they are! In fact, any LFT that maps the interval $[-1, 1]$
to some subinterval $U \subseteq [-1, 1]$ can be interpreted as transforming one number into another number. Quite
surprisingly, and that's the beauty of this representation, we can even prove that we can always produce the next
digit of the result by only looking at a finite - even $O(n)$ in the limit - digits of the operand.

I will leave the details to the [paper][1] but basically, our state is a matrix representing a LFT, and in each step
we either emit a digit $e$ by multiplying from the left with the inverse $M_e^{-1}$ or we shift another digit $d$
from the operand and multiply from the right with $M_d$. A nice fact is that this can all be done using
*integer arithmetic* - though we do need BigInt support, which python offers builtin.

For binary operations, we use a $4 x 2$ matrix representing the (generalized) linear fractional transformation
$L(x, y) = \frac{a xy + c x + e y + g}{b xy + d x + f y + h}$. Emitting a digit is again a multiplication from the
left, but this time we need two operations for shifting digits - either from $x$ or from $y$. But after doing a bit
of arithmetic, you will arrive a formula that looks like two matrix multiplication done alongside/interleaved with
each other. Details, again, in the paper.

# Examples

```
# TODO
```

[1]: https://doi.org/10.1007/3-540-45699-6_5
[2]: https://en.wikipedia.org/wiki/Linear_fractional_transformation
