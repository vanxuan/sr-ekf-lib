import { describe, it, expect } from 'vitest'
import { RingBuf } from '../src/ring-buf'

describe('RingBuf', () => {
  it('pushes and retrieves elements in order', () => {
    const buf = new RingBuf(8)
    for (let i = 0; i < 5; i++) buf.push(i * 10)
    expect(buf.length).toBe(5)
    for (let i = 0; i < 5; i++) expect(buf.get(i)).toBe(i * 10)
  })

  it('overwrites oldest elements when full', () => {
    const buf = new RingBuf(4)
    for (let i = 0; i < 6; i++) buf.push(i)
    expect(buf.length).toBe(4)
    expect(buf.get(0)).toBe(2) // oldest surviving
    expect(buf.get(3)).toBe(5) // newest
  })

  it('shift removes oldest element', () => {
    const buf = new RingBuf(8)
    for (let i = 0; i < 5; i++) buf.push(i)
    expect(buf.shift()).toBe(0)
    expect(buf.length).toBe(4)
    expect(buf.shift()).toBe(1)
    expect(buf.length).toBe(3)
  })

  it('shift returns 0 on empty buffer', () => {
    const buf = new RingBuf(8)
    expect(buf.shift()).toBe(0)
    expect(buf.length).toBe(0)
  })

  it('get works after wrap-around', () => {
    const buf = new RingBuf(4)
    for (let i = 0; i < 6; i++) buf.push(i) // buffer wraps at index 4
    expect(buf.get(0)).toBe(2)
    expect(buf.get(1)).toBe(3)
    expect(buf.get(2)).toBe(4)
    expect(buf.get(3)).toBe(5)
  })

  it('clear resets buffer', () => {
    const buf = new RingBuf(8)
    for (let i = 0; i < 5; i++) buf.push(i)
    buf.clear()
    expect(buf.length).toBe(0)
    expect(buf.shift()).toBe(0)
  })

  it('push updates length up to capacity', () => {
    const buf = new RingBuf(3)
    expect(buf.length).toBe(0)
    buf.push(1); expect(buf.length).toBe(1)
    buf.push(2); expect(buf.length).toBe(2)
    buf.push(3); expect(buf.length).toBe(3)
    buf.push(4); expect(buf.length).toBe(3) // stays at capacity
  })
})
