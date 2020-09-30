// https://www.secg.org/sec2-v2.pdf
// Curve fomula is y^2 = x^3 + ax + b
const CURVE = {
    // Params: a, b
    a: 0n,
    b: 7n,
    // Field over which we'll do calculations
    P: 2n ** 256n - 2n ** 32n - 977n,
    // Subgroup order aka prime_order
    n: 2n ** 256n - 432420386565659656852420866394968145599n,
    // Cofactor
    h: 1n,
    // Base point (x, y) aka generator point
    Gx: 55066263022277343669578718895168534326250603453777594175500187360389116729240n,
    Gy: 32670510020758816978083085130507043184471273380659243275938904335757337482424n,

    // For endomorphism, see below.
    beta: 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501een,
};

class Point {
    static ZERO = new Point(0n, 0n);

    constructor(public x: bigint, public y: bigint) {
    }

    double(): Point {
        const X1 = this.x;
        const Y1 = this.y;
        const lam = mod(3n * X1 ** 2n * invert(2n * Y1, CURVE.P));
        const X3 = mod(lam * lam - 2n * X1);
        const Y3 = mod(lam * (X1 - X3) - Y1);
        return new Point(X3, Y3);
    }

    add(other: Point): Point {
        const [a, b] = [this, other];
        const [X1, Y1, X2, Y2] = [a.x, a.y, b.x, b.y];
        if (X1 == 0n || Y1 == 0n) return b;
        if (X2 == 0n || Y2 == 0n) return a;
        if (X1 == X2 && Y1 == Y2) return this.double();
        if (X1 == X2 && Y1 == -Y2) return Point.ZERO;
        const lam = mod(Y2 - Y1) * invert(X2 - X1, CURVE.P);
        const X3 = mod(lam * lam - X1 - X2);
        const Y3 = mod(lam * (X1 - X3) - Y1);
        return new Point(X3, Y3);
    }

    multiplyDA(n: bigint) {
        let p = Point.ZERO;
        let d: Point = this;
        while (n > 0n) {
            if (n & 1n) p = p.add(d);
            d = d.double();
            n >>= 1n;
        }
        return p;
    }
}

function mod(a: bigint, b: bigint = CURVE.P) {
    const result = a % b;
    return result >= 0 ? result : b + result;
}

// Eucledian GCD
// https://brilliant.org/wiki/extended-euclidean-algorithm/
function egcd(a: bigint, b: bigint) {
    let [x, y, u, v] = [0n, 1n, 1n, 0n];
    while (a !== 0n) {
        let q = b / a;
        let r = b % a;
        let m = x - u * q;
        let n = y - v * q;
        [b, a] = [a, r];
        [x, y] = [u, v];
        [u, v] = [m, n];
    }
    const gcd = b;
    return [gcd, x, y];
}

function invert(number: bigint, modulo: bigint = CURVE.P) {
    if (number == 0n || number <= 0n) {
        throw new Error('invert: expected positive integers');
    }
    let [gcd, x] = egcd(mod(number, modulo), modulo);
    if (gcd !== 1n) {
        throw new Error('invert: does not exist');
    }
    return mod(x, modulo);
}

const G = new Point(CURVE.Gx, CURVE.Gy);

function getPublicKey(privKey: bigint) {
    return G.multiplyDA(privKey);
}

console.log(getPublicKey(140n));