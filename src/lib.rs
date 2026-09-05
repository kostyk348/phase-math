//! # phase-math — обратимая математика для парадигмы PHASE.
//!
//! Математика «внутри фазы» = **биекции на Z/2^64 с аналитически вычисляемой
//! инверсией** (композиция любых = биекция → любая формула точно обратима:
//! `F⁻¹(F(x)) == x`). Float-математика необратима → живёт на границе фазы.
//!
//! - `MulOdd(m)` — умножение на нечётное (биекция; инверсия = mod_inv Newton);
//! - `Add(c)`, `Xor(k)`, `RotL(r)` — базовые биекции;
//! - `arx_round`/`arx_round_inv` — обратимый ARX-раунд;
//! - `Compose` — билдер последовательности с автоматической инверсией;
//! - Q16.16 (`Q`) — фиксированная точка: mul/div/sqrt с **домен-валидацией
//!   на границе** (возвращают `Option` — вызывающий откатывает фазу).

mod biject;
mod q;

pub use biject::*;
pub use q::{Q, q_from_i32, q_mul, q_div, q_sqrt};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compose_roundtrip() {
        let mut c = Compose::new();
        c.push(Op::MulOdd(0x9E3779B97F4A7C15));
        c.push(Op::Add(0xD1B54A32D192ED03));
        c.push(Op::Xor(0xABCDEF0123456789));
        c.push(Op::RotL(19));
        c.push(Op::MulOdd(3));
        let inv = c.inverse();
        let mut rng = 0xF00Du64;
        let mut next = move || {
            rng ^= rng >> 12;
            rng ^= rng << 25;
            rng ^= rng >> 27;
            rng.wrapping_mul(0x2545_F491_4F6C_DD1D)
        };
        for _ in 0..3000 {
            let x = next();
            assert_eq!(inv.apply(c.apply(x)), x);
        }
    }

    #[test]
    fn q_roundtrip_and_domain() {
        // mul: x*y обратим через деление, если y нечётный в Q? Q16.16 mul — не
        // биекция, но домен-валидируется; проверяем точность и границы.
        let a = q_from_i32(3); // 3.0
        let b = q_from_i32(7); // 7.0
        let p = q_mul(a, b).expect("в пределах");
        assert_eq!(p, q_from_i32(21));
        let d = q_div(p, a).expect("ok");
        assert_eq!(d, q_from_i32(7));
        // overflow домена -> None (граница: вызывающий откатывает)
        assert!(q_mul(Q::MAX, Q::MAX).is_none());
        assert!(q_div(a, 0i64.into()).is_none());
        // sqrt точность для полных квадратов
        let s = q_sqrt(q_from_i32(16)).expect("16 -> 4");
        assert_eq!(s, q_from_i32(4));
    }
}
