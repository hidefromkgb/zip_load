#ifndef ZIP_LOAD_H
#define ZIP_LOAD_H

#include <stdint.h>
#include <stdlib.h>

/* [ internal structure, do not use ] */
struct _ZIP_DATA {
    struct {
        uint8_t m_code_size[288]; /* <----------v--- max. Huffman symbols */
        int16_t m_look_up[1 << 10], m_tree[2 * 288]; /*                 | */
    } m_tables[3]; /*          ^-- lookup bits                          | */
    uint8_t m_dict[1 << 15], m_len_codes[137 + 32 + 288]; /* <----------' */
    uint32_t m_num_bits, m_z_adler32, m_final,
             m_check_adler32, m_dist, m_counter, m_num_extra;
    unsigned char *next_in, *next_out;
    unsigned int  avail_in, avail_out;
    size_t m_dist_from_out_buf_start;
    uint64_t m_bit_buf;
};

/* [ internal function, do not use ] */
long _ZIP_Read(struct _ZIP_DATA *data, uint32_t *state, size_t *ilen,
               uint8_t *obgn, uint8_t *onxt, size_t *olen) {
    #define ZIP_CR_RETURN(state_index, result, wait) \
        do {                                         \
            status = result;                         \
            *state = state_index;                    \
            goto common_exit;                        \
            case state_index:;                       \
        } while(wait)

    #define ZIP_GET_BYTE(state_index, c)                              \
        do {                                                          \
            while (pIn_buf_cur >= pIn_buf_end)                        \
                ZIP_CR_RETURN(state_index, (decomp_flags & ZIP_FHMI)? \
                                            ZIP_READ : ZIP_PROG, 0);  \
            c = *pIn_buf_cur++;                                       \
        } while(0)

    #define ZIP_GET_BITS(state_index, b, n)          \
        do {                                         \
            unsigned c;                              \
                                                     \
            while (num_bits < (unsigned)(n)) {       \
                ZIP_GET_BYTE(state_index, c);        \
                bits |= (((ZIP_BITS)c) << num_bits); \
                num_bits += 8;                       \
            }                                        \
            b = bits & ((1 << (n)) - 1);             \
            bits >>= (n);                            \
            num_bits -= (n);                         \
        } while(0)

    #define ZIP_HUFF_DECODE(state_index, sym, pHuff)                         \
        do {                                                                 \
            int temp;                                                        \
            unsigned c, clen;                                                \
                                                                             \
            if (num_bits < 15) {                                             \
                if ((pIn_buf_end - pIn_buf_cur) < 2) {                       \
                    do {                                                     \
                        temp = (pHuff)->m_look_up[bits & (lookup - 1)];      \
                        if (temp >= 0) {                                     \
                            clen = temp >> 9;                                \
                            if ((clen) && (num_bits >= clen))                \
                                break;                                       \
                        }                                                    \
                        else if (num_bits > lookup_bits) {                   \
                            clen = lookup_bits;                              \
                            do temp =                                        \
                            (pHuff)->m_tree[~temp + ((bits >> clen++) & 1)]; \
                            while ((temp < 0) && (num_bits >= (clen + 1)));  \
                            if (temp >= 0)                                   \
                                break;                                       \
                        }                                                    \
                        ZIP_GET_BYTE(state_index, c);                        \
                        bits |= (((ZIP_BITS)c) << num_bits);                 \
                        num_bits += 8;                                       \
                    } while (num_bits < 15);                                 \
                }                                                            \
                else {                                                       \
                    bits |= (((ZIP_BITS)pIn_buf_cur[0]) << (num_bits    ))   \
                         |  (((ZIP_BITS)pIn_buf_cur[1]) << (num_bits + 8));  \
                    pIn_buf_cur += 2;                                        \
                    num_bits += 16;                                          \
                }                                                            \
            }                                                                \
            if ((temp = (pHuff)->m_look_up[bits & (lookup - 1)]) >= 0)       \
                clen = temp >> 9, temp &= 511;                               \
            else {                                                           \
                clen = lookup_bits;                                          \
                do temp = (pHuff)->m_tree[~temp + ((bits >> clen++) & 1)];   \
                while (temp < 0);                                            \
            }                                                                \
            sym = temp;                                                      \
            bits >>= clen;                                                   \
            num_bits -= clen;                                                \
        } while(0)

    #define ZIP_CLEAR_OBJ(obj) memset(&(obj), 0, sizeof(obj))

    #if ZIP_USE_64BIT_BITBUF
    typedef uint64_t ZIP_BITS;
    #define ZIP_BITBUF_SIZE (64)
    #else
    typedef uint32_t ZIP_BITS;
    #define ZIP_BITBUF_SIZE (32)
    #endif

    enum {
        ZIP_PROG = -4, /* cannot progress further */
        ZIP_PARM = -3, /* bad parameter(s) given  */
        ZIP_ADLR = -2, /* ADLER32 mismatch        */
        ZIP_OHNO = -1, /* failure                 */
        ZIP_DONE =  0, /* success                 */
        ZIP_READ =  1, /* needs to read more data */
        ZIP_WRIT =  2, /* has more data to write  */

        ZIP_FPZH =  1, /* flag to parse ZLIB header              */
        ZIP_FHMI =  2, /* flag to show the availability of input */
        ZIP_FNOB =  4, /* flag for non-wrapping output buffer    */
        ZIP_FA32 =  8, /* flag to compute ADLER32                */
    };
    static const int s_length_base[] = {
         3,  4,  5,  6,  7,  8,  9, 10,  11,  13,  15,  17,  19, 23,  27,
        31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258,  0,  0};
    static const int s_length_extra[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, 0, 0};
    static const int s_dist_base[] = {
           1,    2,    3,    4,    5,     7,     9,    13,  17,   25,   33,
          49,   65,   97,  129,  193,   257,   385,   513, 769, 1025, 1537,
        2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,   0,    0};
    static const int s_dist_extra[] = {
        0, 0, 0, 0, 1, 1, 2,  2,  3,  3,  4,  4,  5,  5,  6,
        6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};
    static const uint8_t s_length_dezigzag[] = {
        16, 17, 18,  0, 8,  7, 9,  6, 10,
         5, 11,  4, 12, 3, 13, 2, 14,  1, 15};
    static const int s_min_table_sizes[] = {257, 1, 4};
    const uint32_t decomp_flags = ZIP_FA32 | ZIP_FHMI;
    const int lookup = sizeof( data->m_tables->m_look_up)
                     / sizeof(*data->m_tables->m_look_up);
    const int lookup_bits = 10;

    long status = ZIP_OHNO;
    uint32_t num_bits, dist, counter, num_extra, mtype, zhdr0 = 0, zhdr1 = 0;
    ZIP_BITS bits;
    const uint8_t *pIn_buf_cur = data->next_in,
                  *const pIn_buf_end = data->next_in + *ilen;
    uint8_t *pOut_buf_cur = onxt,
            *const pOut_buf_end = onxt + *olen, rawhdr[4];
    size_t mask = (decomp_flags & ZIP_FNOB)?
                  (size_t)-1 : ((onxt - obgn) + *olen) - 1,
           dist_from_out_buf_start;

    /* Ensure the output buffer's size is a power of 2, unless the  */
    /* output buffer is large enough to hold the entire output file */
    if (((mask + 1) & mask) || (onxt < obgn)) {
        *ilen = *olen = 0;
        return ZIP_PARM;
    }
    num_bits = data->m_num_bits;
    bits = data->m_bit_buf;
    dist = data->m_dist;
    counter = data->m_counter;
    num_extra = data->m_num_extra;
    dist_from_out_buf_start = data->m_dist_from_out_buf_start;
    switch (*state) {
    case 0:
        bits = num_bits = dist = counter = num_extra = zhdr0 = zhdr1 = 0;
        data->m_z_adler32 = data->m_check_adler32 = 1;
        if (decomp_flags & ZIP_FPZH) {
            ZIP_GET_BYTE(1, zhdr0);
            ZIP_GET_BYTE(2, zhdr1);
            counter = (((zhdr0 * 256 + zhdr1) % 31 != 0)
                    ||  (zhdr1 & 32) || ((zhdr0 & 15) != 8));
            if (!(decomp_flags & ZIP_FNOB))
                counter |= (((1U << (8U + (zhdr0 >> 4))) > 32768U)
                        ||  (mask + 1 < (size_t)(1U << (8U + (zhdr0 >> 4)))));
            if (counter)
                ZIP_CR_RETURN(36, ZIP_OHNO, 1);
        }
        do {
            ZIP_GET_BITS(3, data->m_final, 3);
            mtype = data->m_final >> 1;
            if (mtype == 0) {
                ZIP_GET_BITS(5, counter, num_bits & 7);
                for (counter = 0; counter < 4; ++counter)
                    if (num_bits)
                        ZIP_GET_BITS(6, rawhdr[counter], 8);
                    else
                        ZIP_GET_BYTE(7, rawhdr[counter]);
                if ((counter = (rawhdr[0] | (rawhdr[1] << 8)))
                !=  (unsigned)(0xFFFF ^ (rawhdr[2] | (rawhdr[3] << 8))))
                    ZIP_CR_RETURN(39, ZIP_OHNO, 1);
                while ((counter) && (num_bits)) {
                    ZIP_GET_BITS(51, dist, 8);
                    while (pOut_buf_cur >= pOut_buf_end)
                        ZIP_CR_RETURN(52, ZIP_WRIT, 0);
                    *pOut_buf_cur++ = (uint8_t)dist;
                    counter--;
                }
                while (counter) {
                    size_t m, n;

                    while (pOut_buf_cur >= pOut_buf_end)
                        ZIP_CR_RETURN(9, ZIP_WRIT, 0);
                    while (pIn_buf_cur >= pIn_buf_end)
                        ZIP_CR_RETURN(38, (decomp_flags & ZIP_FHMI)? ZIP_READ : ZIP_PROG, 0);
                    m = (size_t)(pOut_buf_end - pOut_buf_cur);
                    n = (size_t)(pIn_buf_end - pIn_buf_cur);
                    m = (m < n)? m : n;
                    n = (m < counter)? m : counter;
                    memcpy(pOut_buf_cur, pIn_buf_cur, n);
                    pIn_buf_cur += n;
                    pOut_buf_cur += n;
                    counter -= (unsigned)n;
                }
            } else if (mtype == 3) {
                ZIP_CR_RETURN(10, ZIP_OHNO, 1);
            } else {
                uint32_t table_sizes[sizeof(data->m_tables) / sizeof(*data->m_tables)];

                if (mtype == 1) {
                    uint8_t *p = data->m_tables[0].m_code_size;
                    unsigned i;
                    table_sizes[0] = 288;
                    table_sizes[1] = 32;
                    memset(data->m_tables[1].m_code_size, 5, 32);
                    for (i = 0; i <= 143; ++i)
                        *p++ = 8;
                    for (; i <= 255; ++i)
                        *p++ = 9;
                    for (; i <= 279; ++i)
                        *p++ = 7;
                    for (; i <= 287; ++i)
                        *p++ = 8;
                } else {
                    for (counter = 0; counter < 3; counter++) {
                        ZIP_GET_BITS(11, table_sizes[counter], "\05\05\04"[counter]);
                        table_sizes[counter] += s_min_table_sizes[counter];
                    }
                    ZIP_CLEAR_OBJ(data->m_tables[2].m_code_size);
                    for (counter = 0; counter < table_sizes[2]; counter++) {
                        unsigned s;

                        ZIP_GET_BITS(14, s, 3);
                        data->m_tables[2].m_code_size[s_length_dezigzag[counter]] = (uint8_t)s;
                    }
                    table_sizes[2] = 19;
                }
                for (; (int)mtype >= 0; mtype--) {
                    int tree_next, tree_cur;
                    unsigned i, j, used_syms, total, sym_index, next_code[17], total_syms[16];

                    ZIP_CLEAR_OBJ(total_syms);
                    ZIP_CLEAR_OBJ(data->m_tables[mtype].m_look_up);
                    ZIP_CLEAR_OBJ(data->m_tables[mtype].m_tree);
                    for (i = 0; i < table_sizes[mtype]; ++i)
                        total_syms[data->m_tables[mtype].m_code_size[i]]++;
                    used_syms = 0, total = 0;
                    next_code[0] = next_code[1] = 0;
                    for (i = 1; i <= 15; ++i) {
                        used_syms += total_syms[i];
                        next_code[i + 1] = (total = ((total + total_syms[i]) << 1));
                    }
                    if ((65536 != total) && (used_syms > 1))
                        ZIP_CR_RETURN(35, ZIP_OHNO, 1);
                    for (tree_next = -1, sym_index = 0; sym_index < table_sizes[mtype]; ++sym_index) {
                        unsigned rev_code = 0, l, cur_code, code_size = data->m_tables[mtype].m_code_size[sym_index];
                        if (!code_size)
                            continue;
                        cur_code = next_code[code_size]++;
                        for (l = code_size; l > 0; l--, cur_code >>= 1)
                            rev_code = (rev_code << 1) | (cur_code & 1);
                        if (code_size <= lookup_bits) {
                            int16_t k = (int16_t)((code_size << 9) | sym_index);
                            while (rev_code < lookup) {
                                data->m_tables[mtype].m_look_up[rev_code] = k;
                                rev_code += (1 << code_size);
                            }
                            continue;
                        }
                        if (0 == (tree_cur = data->m_tables[mtype].m_look_up[rev_code & (lookup - 1)])) {
                            data->m_tables[mtype].m_look_up[rev_code & (lookup - 1)] = (int16_t)tree_next;
                            tree_cur = tree_next;
                            tree_next -= 2;
                        }
                        rev_code >>= (lookup_bits - 1);
                        for (j = code_size; j > (lookup_bits + 1); j--) {
                            tree_cur -= ((rev_code >>= 1) & 1);
                            if (!data->m_tables[mtype].m_tree[-tree_cur - 1]) {
                                data->m_tables[mtype].m_tree[-tree_cur - 1] = (int16_t)tree_next;
                                tree_cur = tree_next;
                                tree_next -= 2;
                            } else
                                tree_cur = data->m_tables[mtype].m_tree[-tree_cur - 1];
                        }
                        tree_cur -= ((rev_code >>= 1) & 1);
                        data->m_tables[mtype].m_tree[-tree_cur - 1] = (int16_t)sym_index;
                    }
                    if (mtype == 2) {
                        for (counter = 0; counter < (table_sizes[0] + table_sizes[1]);) {
                            unsigned s;

                            ZIP_HUFF_DECODE(16, dist, &data->m_tables[2]);
                            if (dist < 16) {
                                data->m_len_codes[counter++] = (uint8_t)dist;
                                continue;
                            }
                            if ((dist == 16) && (!counter))
                                ZIP_CR_RETURN(17, ZIP_OHNO, 1);
                            num_extra = "\02\03\07"[dist - 16];
                            ZIP_GET_BITS(18, s, num_extra);
                            s += "\03\03\013"[dist - 16];
                            memset(data->m_len_codes + counter, (dist == 16)? data->m_len_codes[counter - 1] : 0, s);
                            counter += s;
                        }
                        if ((table_sizes[0] + table_sizes[1]) != counter)
                            ZIP_CR_RETURN(21, ZIP_OHNO, 1);
                        memcpy(data->m_tables[0].m_code_size, data->m_len_codes, table_sizes[0]);
                        memcpy(data->m_tables[1].m_code_size, data->m_len_codes + table_sizes[0], table_sizes[1]);
                    }
                }
                while (!0) {
                    uint8_t *pSrc;
                    while (!0) {
                        if (((pIn_buf_end - pIn_buf_cur) < 4) || ((pOut_buf_end - pOut_buf_cur) < 2)) {
                            ZIP_HUFF_DECODE(23, counter, &data->m_tables[0]);
                            if (counter >= 256)
                                break;
                            while (pOut_buf_cur >= pOut_buf_end)
                                ZIP_CR_RETURN(24, ZIP_WRIT, 0);
                            *pOut_buf_cur++ = (uint8_t)counter;
                        } else {
                            int sym2;
                            unsigned clen; /* code length */
#if ZIP_USE_64BIT_BITBUF
                            if (num_bits < 30) {
                                bits |= (((ZIP_BITS)MZ_READ_LE32(pIn_buf_cur)) << num_bits);
                                pIn_buf_cur += 4;
                                num_bits += 32;
                            }
#else
                            if (num_bits < 15) {
                                bits |= (((ZIP_BITS)*((const uint16_t *)(pIn_buf_cur))) << num_bits);
                                pIn_buf_cur += 2;
                                num_bits += 16;
                            }
#endif
                            if ((sym2 = data->m_tables[0].m_look_up[bits & (lookup - 1)]) >= 0)
                                clen = sym2 >> 9;
                            else {
                                clen = lookup_bits;
                                do sym2 = data->m_tables[0].m_tree[~sym2 + ((bits >> clen++) & 1)];
                                while (sym2 < 0);
                            }
                            counter = sym2;
                            bits >>= clen;
                            num_bits -= clen;
                            if (counter & 256)
                                break;

#if !ZIP_USE_64BIT_BITBUF
                            if (num_bits < 15) {
                                bits |= (((ZIP_BITS)*((const uint16_t *)(pIn_buf_cur))) << num_bits);
                                pIn_buf_cur += 2;
                                num_bits += 16;
                            }
#endif
                            if ((sym2 = data->m_tables[0].m_look_up[bits & (lookup - 1)]) >= 0)
                                clen = sym2 >> 9;
                            else {
                                clen = lookup_bits;
                                do sym2 = data->m_tables[0].m_tree[~sym2 + ((bits >> clen++) & 1)];
                                while (sym2 < 0);
                            }
                            bits >>= clen;
                            num_bits -= clen;

                            pOut_buf_cur[0] = (uint8_t)counter;
                            if (sym2 & 256) {
                                pOut_buf_cur++;
                                counter = sym2;
                                break;
                            }
                            pOut_buf_cur[1] = (uint8_t)sym2;
                            pOut_buf_cur += 2;
                        }
                    }
                    if ((counter &= 511) == 256)
                        break;

                    num_extra = s_length_extra[counter - 257];
                    counter = s_length_base[counter - 257];
                    if (num_extra) {
                        unsigned extra_bits;

                        ZIP_GET_BITS(25, extra_bits, num_extra);
                        counter += extra_bits;
                    }

                    ZIP_HUFF_DECODE(26, dist, &data->m_tables[1]);
                    num_extra = s_dist_extra[dist];
                    dist = s_dist_base[dist];
                    if (num_extra) {
                        unsigned extra_bits;

                        ZIP_GET_BITS(27, extra_bits, num_extra);
                        dist += extra_bits;
                    }

                    dist_from_out_buf_start = pOut_buf_cur - obgn;
                    if ((dist > dist_from_out_buf_start) && (decomp_flags & ZIP_FNOB))
                        ZIP_CR_RETURN(37, ZIP_OHNO, 1);

                    pSrc = obgn + ((dist_from_out_buf_start - dist) & mask);

                    if (((pOut_buf_cur > pSrc)? pOut_buf_cur : pSrc) > pOut_buf_end - counter) {
                        while (counter--) {
                            while (pOut_buf_cur >= pOut_buf_end)
                                ZIP_CR_RETURN(53, ZIP_WRIT, 0);
                            *pOut_buf_cur++ = obgn[(dist_from_out_buf_start++ - dist) & mask];
                        }
                        continue;
                    }
#if ZIP_USE_UNALIGNED_LOADS_AND_STORES
                    else if ((counter >= 9) && (counter <= dist)) {
                        const uint8_t *pSrc_end = pSrc + (counter & ~7);
                        do {
#ifdef ZIP_UNALIGNED_USE_MEMCPY
                            memcpy(pOut_buf_cur, pSrc, sizeof(uint32_t)*2);
#else
                            ((uint32_t *)pOut_buf_cur)[0] = ((const uint32_t *)pSrc)[0];
                            ((uint32_t *)pOut_buf_cur)[1] = ((const uint32_t *)pSrc)[1];
#endif
                            pOut_buf_cur += 8;
                        } while ((pSrc += 8) < pSrc_end);
                        if ((counter &= 7) < 3) {
                            if (counter) {
                                pOut_buf_cur[0] = pSrc[0];
                                if (counter > 1)
                                    pOut_buf_cur[1] = pSrc[1];
                                pOut_buf_cur += counter;
                            }
                            continue;
                        }
                    }
#endif
                    while(counter>2) {
                        pOut_buf_cur[0] = pSrc[0];
                        pOut_buf_cur[1] = pSrc[1];
                        pOut_buf_cur[2] = pSrc[2];
                        pOut_buf_cur += 3;
                        pSrc += 3;
                        counter -= 3;
                    }
                    if (counter > 0) {
                        pOut_buf_cur[0] = pSrc[0];
                        if (counter > 1)
                            pOut_buf_cur[1] = pSrc[1];
                        pOut_buf_cur += counter;
                    }
                }
            }
        } while (!(data->m_final & 1));

        /* Ensure byte alignment and put back any bytes from the bitbuf if */
        /* we have looked ahead too far on gzip, or other Deflate streams  */
        /* followed by arbitrary data. I'm being super conservative here.  */
        /* A number of simplifications can be made to the byte alignment   */
        /* part, and the Adler32 check shouldn't ever need to worry about  */
        /* reading from the bit buffer now. */
        uint32_t null;

        ZIP_GET_BITS(32, null, num_bits & 7);
        for (null = 8; (pIn_buf_cur > data->next_in) && (num_bits >= null);
             num_bits -= null, --pIn_buf_cur);
        bits &= (ZIP_BITS)((((uint64_t)1) << num_bits) - (uint64_t)1);

        if (decomp_flags & ZIP_FPZH) {
            for (counter = 0; counter < 4; ++counter) {
                unsigned s;

                if (num_bits)
                    ZIP_GET_BITS(41, s, 8);
                else
                    ZIP_GET_BYTE(42, s);
                data->m_z_adler32 = (data->m_z_adler32 << 8) | s;
            }
        }
        ZIP_CR_RETURN(34, ZIP_DONE, 1);
    }

common_exit:
    /* As long as we aren't telling the caller that we NEED more input to   */
    /* make further progress: put back any bytes from the bit buffer if we  */
    /* have looked ahead too far on gzip, or other Deflate streams followed */
    /* by arbitrary data. We need to be very careful here though, to NOT    */
    /* push back any bytes we definitely know we need to make further       */
    /* progress, or we will end up locking the caller in an infinite loop.  */
    if ((status != ZIP_READ) && (status != ZIP_PROG)) {
        while ((pIn_buf_cur > data->next_in) && (num_bits >= 8)) {
            --pIn_buf_cur;
            num_bits -= 8;
        }
    }
    data->m_num_bits = num_bits;
    data->m_bit_buf = bits & (ZIP_BITS)((1LLU << num_bits) - 1LLU);
    data->m_dist = dist;
    data->m_counter = counter;
    data->m_num_extra = num_extra;
    data->m_dist_from_out_buf_start = dist_from_out_buf_start;
    *ilen = pIn_buf_cur - data->next_in;
    *olen = pOut_buf_cur - onxt;
    if ((decomp_flags & (ZIP_FPZH | ZIP_FA32)) && (status >= 0)) {
        const uint8_t *ptr = onxt;
        size_t buf_len = *olen;
        uint32_t i, s1 = (uint16_t)data->m_check_adler32,
                    s2 = data->m_check_adler32 >> 16;
        size_t block_len = buf_len % 5552;
        while (buf_len) {
            for (i = 0; i + 7 < block_len; i += 8, ptr += 8) {
                s1 += ptr[0], s2 += s1; s1 += ptr[1], s2 += s1;
                s1 += ptr[2], s2 += s1; s1 += ptr[3], s2 += s1;
                s1 += ptr[4], s2 += s1; s1 += ptr[5], s2 += s1;
                s1 += ptr[6], s2 += s1; s1 += ptr[7], s2 += s1;
            }
            for (; i < block_len; ++i)
                s1 += *ptr++, s2 += s1;
            s1 %= 65521U, s2 %= 65521U;
            buf_len -= block_len;
            block_len = 5552;
        }
        data->m_check_adler32 = (s2 << 16) + s1;
        if ((status == ZIP_DONE) && (decomp_flags & ZIP_FPZH)
        &&  (data->m_check_adler32 != data->m_z_adler32))
            status = ZIP_ADLR;
    }
    return status;
    #undef ZIP_CLEAR_OBJ
    #undef ZIP_HUFF_DECODE
    #undef ZIP_GET_BITS
    #undef ZIP_GET_BYTE
    #undef ZIP_CR_RETURN
}

void ZIP_Load(char *file, long size, void *user,
              void (*save)(char*, char*, long, void*)) {
#pragma pack(push, 1)
    struct {
        uint32_t head; /* header                  */
        uint16_t xver, /* compressor version      */
                 flgs, /* general-purpose flags   */
                 func, /* compression function    */
                 time, /* last modification time  */
                 date; /* last modification date  */
        uint32_t xcrc, /* CRC32 sum               */
                 szcp, /* size, compressed        */
                 szun; /* size, uncompressed      */
        uint16_t szfn, /* size of the file name   */
                 szxf; /* size of the extra field */
    } *zhdr = (void*)file;
#pragma pack(pop)
    uint8_t *retn, *rtmp, *halt = (uint8_t*)file + size;
    struct _ZIP_DATA *data = (void*)malloc(sizeof(*data));
    long stat = 0;

    for (; (uint8_t*)zhdr < halt; zhdr = (void*)
        ((uint8_t*)(zhdr + 1) + zhdr->szfn + zhdr->szxf + zhdr->szcp)) {
        file = (zhdr->szfn && (zhdr->head == 0x04034B50))?
                memcpy(calloc(zhdr->szfn + 1, 1), zhdr + 1, zhdr->szfn) : 0;
        if (!zhdr->szcp && file)
            save(file, 0, 0, user);
        else if (zhdr->szcp && file) {
            rtmp = retn = (uint8_t*)(zhdr + 1) + zhdr->szfn + zhdr->szxf;
            if (zhdr->func == 0x08) {
                unsigned n, dict_ofs = 0, dict = 0;
                size_t i_sz, out_bytes, orig_avail_in;
                uint32_t state = 0;

                data->next_in = rtmp;
                data->avail_in = orig_avail_in = zhdr->szcp;
                data->next_out = retn = calloc(1, zhdr->szun);
                data->avail_out = zhdr->szun;
                while (!0) {
                    i_sz = data->avail_in;
                    out_bytes = sizeof(data->m_dict) - dict_ofs;
                    stat = _ZIP_Read(data, &state, &i_sz, data->m_dict,
                                     data->m_dict + dict_ofs, &out_bytes);
                    data->next_in += (unsigned)i_sz;
                    data->avail_in -= (unsigned)i_sz;
                    dict = (unsigned)out_bytes; /* <-- available dictionary */
                    n = (dict < data->avail_out)? dict : data->avail_out;
                    memcpy(data->next_out, data->m_dict + dict_ofs, n);
                    data->next_out += n;
                    data->avail_out -= n;
                    dict_ofs = (dict_ofs + n) & (sizeof(data->m_dict) - 1);
                    dict -= n;
                    if (stat < 0) {
                        stat = -1; /* Stream is corrupt (maybe some data was */
                        break;     /* left uncompressed in the output dict?) */
                    }
                    else if ((stat == 1) && !orig_avail_in) {
                        stat = -2; /* Cannot make further progress */
                        break;     /* without supplying more input */
                    }
                    else if (!data->avail_in  || !stat
                         ||  !data->avail_out ||  dict) {
                        stat = (!stat && !dict)? 1 : 0;
                        break; /* 1 = stream end; 0 = success */
                    }
                }
            }
            if (stat >= 0)
                save(file, (char*)retn, zhdr->szun, user);
            if (zhdr->func == 0x08)
                free(retn);
        }
        free(file);
    }
    free(data);
}

#endif /* ZIP_LOAD_H */
