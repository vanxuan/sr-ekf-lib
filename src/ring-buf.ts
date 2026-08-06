export class RingBuf {
  private buf: Float64Array;
  private head = 0;
  length = 0;
  private readonly mask: number;

  constructor(capacity: number) {
    this.buf = new Float64Array(capacity);
    this.mask = capacity - 1;
  }

  push(v: number): void {
    this.buf[this.head] = v;
    this.head = (this.head + 1) & this.mask;
    if (this.length < this.buf.length) this.length++;
  }

  shift(): number {
    if (this.length === 0) return 0;
    const idx = (this.head - this.length + this.buf.length) & this.mask;
    this.length--;
    return this.buf[idx];
  }

  get(i: number): number {
    const idx = (this.head - this.length + i + this.buf.length) & this.mask;
    return this.buf[idx];
  }

  clear(): void {
    this.head = 0;
    this.length = 0;
  }
}
