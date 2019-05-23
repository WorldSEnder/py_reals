# Computing with real numbers (in python)

The implementation is based on [Computing with Real Numbers][1]. Here's a quick breakdown for
the interested reader:

As our primitive, we model only numbers in the interval <img src="/svgs/43ca5ad9e1f094a31392f860ef481e5c.svg" align=middle width=45.66218414999998pt height=24.65753399999998pt/>. Other numbers will be represented
by scaling appropriately. We will work in base <img src="/svgs/76c5792347bb90ef71cfbace628572cf.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/>, although this can be tweaked in the implementation.
Any power of two should work fine out of the box. A number for us will be an infinite stream of digits.

Consider an ordinary digit stream of base <img src="/svgs/76c5792347bb90ef71cfbace628572cf.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/>, meaning you either get a <img src="/svgs/29632a9bf827ce0200454dd32fc3be82.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> or a <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/>. The normal arithmetic
tells you that seeing a <img src="/svgs/29632a9bf827ce0200454dd32fc3be82.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> as the first digit means the rest of the number must be in the interval <img src="/svgs/616b42c56c829ade05e78d7e21bbd101.svg" align=middle width=47.48867639999999pt height=24.65753399999998pt/>,
seeing a <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> similarly restricts to <img src="/svgs/c3b36e706022f485d4e484b9889bf969.svg" align=middle width=47.48867639999999pt height=24.65753399999998pt/>. The problem with this approach is obvious once you try
to add two such numbers. When one of the numbers is a stream of <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/>s and the other one is a stream of <img src="/svgs/29632a9bf827ce0200454dd32fc3be82.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> you
still can not possibly say if the result must begin with a <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/>! Any <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> in the second stream would lead to
a carry and the first digit of the result would be a <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/>. If there is no <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> in the second stream, there would be
no carry, and the first digit of the result would be a <img src="/svgs/29632a9bf827ce0200454dd32fc3be82.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/>.

As a remedy, we will allow *negative* digits, also, which we write as <img src="/svgs/6c1a1040199f9f849926e89ad86ae1ec.svg" align=middle width=8.219209349999991pt height=25.698634500000015pt/>.
Formally, we will interpret a stream of digits as a so called [Linear fractional transformation][2], say <img src="/svgs/bc8afb5802dbe3a960f6436ffa035ed4.svg" align=middle width=33.36766454999999pt height=24.65753399999998pt/>,
and after seeing a finite number of digits, the real number of the whole stream is somewhere in the interval
<img src="/svgs/3e1572360e205694afb6ed150ea65c6f.svg" align=middle width=69.63485924999999pt height=24.65753399999998pt/>. To be more precise, we will interpret a <img src="/svgs/6c1a1040199f9f849926e89ad86ae1ec.svg" align=middle width=8.219209349999991pt height=25.698634500000015pt/> as the transformation
<img src="/svgs/334ca11295fd8292f5d77245cadfe4e1.svg" align=middle width=93.67356405pt height=27.77565449999998pt/>, <img src="/svgs/29632a9bf827ce0200454dd32fc3be82.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> is interpreted as <img src="/svgs/63386a56e9b7ec91d93476d499ac9fa0.svg" align=middle width=76.84699604999999pt height=24.65753399999998pt/> and <img src="/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align=middle width=8.219209349999991pt height=21.18721440000001pt/> as <img src="/svgs/79f73098b6ce10ef15968e9fd1fea95b.svg" align=middle width=93.49091729999999pt height=27.77565449999998pt/>.
Seeing multiple digits after another <img src="/svgs/465a8cd3ee4c5ce3061192c40af6d51a.svg" align=middle width=16.21010159999999pt height=22.831056599999986pt/> is simply the composition of those functions <img src="/svgs/53488d47e9e0293853a58024ba7f1d4a.svg" align=middle width=171.05712194999998pt height=24.65753399999998pt/>.

## Calculation of compositions of linear fractional transformations

Before we go to calculating anything, there is a nice trick for calculating with LFTs. Namely, if you have
two LFTs <img src="/svgs/a9fb980bcf6736e1a873e92f843b3406.svg" align=middle width=88.85641379999998pt height=28.92634470000001pt/> and <img src="/svgs/946943a34d16d2a06800f2b30b33d440.svg" align=middle width=93.80870564999998pt height=27.19121129999998pt/>, then their composition can be calculated
easily. Suggestively reading each LFT as a matrix, their composition is just the normal matrix product!

<p><img src="/svgs/4557afcd1ff1bdfa12f19524b215141e.svg" align=middle width=416.69056605pt height=137.07763245pt/></p>

## A simple example of a single number

Let's do an example to get warmed up: Say, we observe the stream of digits <img src="/svgs/33f8b266959943cc51a80d090e114655.svg" align=middle width=48.40183589999999pt height=25.698634500000015pt/>, then the number will
be contained somewhere in the following interval:

<p><img src="/svgs/beb1d1765adfc79b4ea6b61393487603.svg" align=middle width=246.01606259999997pt height=195.89433269999998pt/></p>

We conclude that the number starting with <img src="/svgs/33f8b266959943cc51a80d090e114655.svg" align=middle width=48.40183589999999pt height=25.698634500000015pt/> must be somewhere in the interval <img src="/svgs/3e35ff747f4b85f4e398f985fd756aa0.svg" align=middle width=66.66682604999998pt height=24.65753399999998pt/>.

## Unary and binary operations

Since we now know how to interpret a single number, we can see if we can build up some operations on them. It would
be convenient, if those were also LFTs and indeed, they are! In fact, any LFT that maps the interval <img src="/svgs/43ca5ad9e1f094a31392f860ef481e5c.svg" align=middle width=45.66218414999998pt height=24.65753399999998pt/>
to some subinterval <img src="/svgs/9b0b7334753ccadbe1f9494f91e457ba.svg" align=middle width=80.59575314999998pt height=24.65753399999998pt/> can be interpreted as transforming one number into another number. Quite
surprisingly, and that's the beauty of this representation, we can even prove that we can always produce the next
digit of the result by only looking at a finite - even <img src="/svgs/1f08ccc9cd7309ba1e756c3d9345ad9f.svg" align=middle width=35.64773519999999pt height=24.65753399999998pt/> in the limit - digits of the operand.

I will leave the details to the [paper][1] but basically, our state is a matrix representing a LFT, and in each step
we either emit a digit <img src="/svgs/8cd34385ed61aca950a6b06d09fb50ac.svg" align=middle width=7.654137149999991pt height=14.15524440000002pt/> by multiplying from the left with the inverse <img src="/svgs/b36f150c18d3afade293a22dff116fba.svg" align=middle width=34.56628394999999pt height=26.76175259999998pt/> or we shift another digit <img src="/svgs/2103f85b8b1477f430fc407cad462224.svg" align=middle width=8.55596444999999pt height=22.831056599999986pt/>
from the operand and multiply from the right with <img src="/svgs/e4a4ace2d727b81133d571213c643be9.svg" align=middle width=22.79059199999999pt height=22.465723500000017pt/>. A nice fact is that this can all be done using
*integer arithmetic* - though we do need BigInt support, which python offers builtin.

For binary operations, we use a <img src="/svgs/68dfa78cbdac03942b1b75f702cf488b.svg" align=middle width=25.833406499999988pt height=21.18721440000001pt/> matrix representing the (generalized) linear fractional transformation
<img src="/svgs/64b0b9ff48758b60067ed39eb429e54d.svg" align=middle width=160.57486335pt height=28.913154600000016pt/>. Emitting a digit is again a multiplication from the
left, but this time we need two operations for shifting digits - either from <img src="/svgs/332cc365a4987aacce0ead01b8bdcc0b.svg" align=middle width=9.39498779999999pt height=14.15524440000002pt/> or from <img src="/svgs/deceeaf6940a8c7a5a02373728002b0f.svg" align=middle width=8.649225749999989pt height=14.15524440000002pt/>. But after doing a bit
of arithmetic, you will arrive a formula that looks like two matrix multiplication done alongside/interleaved with
each other. Details, again, in the paper.

# Examples

```
# TODO
```

[1]: https://doi.org/10.1007/3-540-45699-6_5
[2]: https://en.wikipedia.org/wiki/Linear_fractional_transformation
