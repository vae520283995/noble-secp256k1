/*! noble-secp256k1 - MIT License (c) Paul Miller (paulmillr.com) */
'use strict';
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

let precomputes: Point[];
let precomputesJ: JacobianPoint[];

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
        if (X1 === 0n || Y1 === 0n) return b;
        if (X2 === 0n || Y2 === 0n) return a;
        if (X1 === X2 && Y1 === Y2) return this.double();
        if (X1 === X2 && Y1 === -Y2) return Point.ZERO;
        const lam = mod((Y2 - Y1) * invert(X2 - X1, CURVE.P));
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

    getPrecomputes(): Point[] {
        if (precomputes) return precomputes;
        precomputes = [];
        let dbl: Point = this;
        for (let i = 0; i < 256; i++) {
            precomputes.push(dbl);
            dbl = dbl.double();
        }
        return precomputes;
    }

    multiplyCT(n: bigint) {
        let dbl = new Point(this.x, this.y);
        let p = Point.ZERO;
        let f = Point.ZERO;
        let dbls = this.getPrecomputes();
        for (let i = 0; dbl = dbls[i], i < 256; i++) {
            if ((n & 1n) === 1n) p = p.add(dbl);
            else f = f.add(dbl);
            n >>= 1n;
        }
        return p;
    }

    getPrecomputesJ(): JacobianPoint[] {
        if (precomputesJ) return precomputesJ;
        precomputesJ = [];
        let dbl: JacobianPoint = JacobianPoint.fromAffine(this);
        for (let i = 0; i < 256; i++) {
            precomputesJ.push(dbl);
            dbl = dbl.double();
        }
        return precomputesJ;
    }

    multiplyPreCTJ(n: bigint) {
        const precomputesJ = this.getPrecomputesJ();
        let dbl = JacobianPoint.fromAffine(this);
        let p = JacobianPoint.ZERO;
        let f = JacobianPoint.ZERO;
        for (let i = 0; dbl = precomputesJ[i], i < 256; i++) {
            if ((n & 1n)  === 1n) p = p.add(dbl);
            else f = f.add(dbl);
            n >>= 1n;
        }
        return p.toAffine();
    }
}

class JacobianPoint {
    constructor(public x: bigint, public y: bigint, public z: bigint) {}

    static BASE = new JacobianPoint(CURVE.Gx, CURVE.Gy, 1n);
    static ZERO = new JacobianPoint(0n, 1n, 0n);

    static fromAffine(p: Point): JacobianPoint {
        if (!(p instanceof Point)) {
            throw new TypeError('JacobianPoint#fromAffine: expected Point');
        }
        return new JacobianPoint(p.x, p.y, 1n);
    }

    toAffine(invZ: bigint = invert(this.z)): Point {
        const invZ2 = invZ ** 2n;
        const x = mod(this.x * invZ2);
        const y = mod(this.y * invZ2 * invZ);
        return new Point(x, y);
    }

    // Fast algo for doubling 2 Jacobian Points when curve's a=0.
    // Note: cannot be reused for other curves when a != 0.
    // From: http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
    // Cost: 2M + 5S + 6add + 3*2 + 1*3 + 1*8.
    double(): JacobianPoint {
        const X1 = this.x;
        const Y1 = this.y;
        const Z1 = this.z;
        const A = X1 ** 2n;
        const B = Y1 ** 2n;
        const C = B ** 2n;
        const D = 2n * ((X1 + B) ** 2n - A - C);
        const E = 3n * A;
        const F = E ** 2n;
        const X3 = mod(F - 2n * D);
        const Y3 = mod(E * (D - X3) - 8n * C);
        const Z3 = mod(2n * Y1 * Z1);
        return new JacobianPoint(X3, Y3, Z3);
    }

    // Fast algo for adding 2 Jacobian Points when curve's a=0.
    // Note: cannot be reused for other curves when a != 0.
    // http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-1998-cmo-2
    // Cost: 12M + 4S + 6add + 1*2.
    // Note: 2007 Bernstein-Lange (11M + 5S + 9add + 4*2) is actually *slower*. No idea why.
    add(other: JacobianPoint): JacobianPoint {
        if (!(other instanceof JacobianPoint)) {
            throw new TypeError('JacobianPoint#add: expected JacobianPoint');
        }
        const X1 = this.x;
        const Y1 = this.y;
        const Z1 = this.z;
        const X2 = other.x;
        const Y2 = other.y;
        const Z2 = other.z;
        if (X2 === 0n || Y2 === 0n) return this;
        if (X1 === 0n || Y1 === 0n) return other;
        const Z1Z1 = Z1 ** 2n;
        const Z2Z2 = Z2 ** 2n;
        const U1 = X1 * Z2Z2;
        const U2 = X2 * Z1Z1;
        const S1 = Y1 * Z2 * Z2Z2;
        const S2 = Y2 * Z1 * Z1Z1;
        const H = mod(U2 - U1);
        const r = mod(S2 - S1);
        // H = 0 meaning it's the same point.
        if (H === 0n) {
            if (r === 0n) {
                return this.double();
            } else {
                return JacobianPoint.ZERO;
            }
        }
        const HH = mod(H ** 2n);
        const HHH = mod(H * HH);
        const V = U1 * HH;
        const X3 = mod(r ** 2n - HHH - 2n * V);
        const Y3 = mod(r * (V - X3) - S1 * HHH);
        const Z3 = mod(Z1 * Z2 * H);
        return new JacobianPoint(X3, Y3, Z3);
    }
}

function mod(a: bigint, b: bigint = CURVE.P): bigint {
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
    if (number === 0n || modulo <= 0n) {
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

console.log(getPublicKey(2n ** 255n - 19n));
console.log(G.multiplyPreCTJ(2n ** 255n - 19n));

const {run, mark, logMem} = require('micro-bmark');

// run([4, 8, 16], async (windowSize) => {
run(async () => {
    const samples = 1;
    //console.log(`-------\nBenchmarking window=${windowSize} samples=${samples}...`);

    logMem();
    console.log();

    // await mark('getPublicKey 256 bit', samples, () => {
    //     G.multiplyCT(2n ** 255n - 42n);
    // });
    //
    // await mark('getPublicKey 1 bit', samples * 10, () => {
    //     G.multiplyCT(2n);
    // });
    //
    // await mark('getPublicKey 256 bit', samples * 10, () => {
    //     G.multiplyCT(2n ** 255n - 42n);
    // });

    await mark('multiplyPreCTJ 256 bit', samples, () => {
        G.multiplyPreCTJ(2n ** 255n - 42n);
    });

    await mark('multiplyPreCTJ 1 bit', samples * 10, () => {
        G.multiplyPreCTJ(2n);
    });

    await mark('multiplyPreCTJ 256 bit', samples * 10, () => {
        G.multiplyPreCTJ(2n ** 255n - 42n);
    });

    console.log();
    logMem();
});

