/* XSPICE code model for the Controlled PWM Oscillator.
 * This is a complete redesign of the original version,
 * according to the d_osc model provided by G. Atkinson
 */

#include <stdlib.h>

// Controls for timing of next scheduled call. */

#define FACTOR1 0.75
#define FACTOR2 0.8

/* PWL table entry. */

struct pwl {
    double ctl, dc;
};

/* Called at end to free memory. */

static void cm_d_pwm_callback(ARGS, Mif_Callback_Reason_t reason)
{
    if (reason == MIF_CB_DESTROY) {
        struct pwl *table = STATIC_VAR(locdata);

        if (table)
            free(table);
        STATIC_VAR(locdata) = NULL;
    }
}

/* Get the current duty cycle. */

static double get_dc(double ctl, struct pwl *table, int csize)
{
    double d;
    int    i;

    for (i = 0; i < csize; ++i) {
         if (table[i].ctl > ctl)
             break;
    }

    /* Interpolation outside input range continues slope. */

    if (i > 0) {
       if (i == csize)
           i -= 2;
       else
           i--;
    }
    d = table[i].dc +
            (ctl - table[i].ctl) * ((table[i + 1].dc - table[i].dc) /
                                    (table[i + 1].ctl - table[i].ctl));

    /* limit duty cycle d to 0.01 <= d <= 0.99 */
    if (d > 0.99)
        d = 0.99;
    else if (d < 0.01)
        d = 0.01;

    return d;
}

/* The state data. */

struct state {
    double          last_time; // Time of last output change.
    Digital_State_t last;      // Last value output.
};

/* The code-model function. */

void cm_d_pwm(ARGS)
{
    struct pwl   *table;
    struct state *state;
    double        ctl, delta, when, ddc;
    int           csize, i;

    CALLBACK = cm_d_pwm_callback;

    csize = PARAM_SIZE(cntl_array);
    delta = 1.0 / PARAM(frequency);
    
    if (INIT) {

        /* Validate PWL table. */

        for (i = 0; i < csize - 1; ++i) {
            if (PARAM(cntl_array[i]) >= PARAM(cntl_array[i + 1]))
                break;
        }

        if (i < csize - 1 || csize != PARAM_SIZE(dc_array)) {
            cm_message_send("Badly-formed control table");
            STATIC_VAR(locdata) = NULL;
            return;
        }

        /* Allocate PWL table. */

        table = malloc(csize * sizeof (struct pwl));
        STATIC_VAR(locdata) = table;
        if (!table)
            return;

        for (i = 0; i < csize; ++i) {
            table[i].ctl = PARAM(cntl_array[i]);
            table[i].dc = PARAM(dc_array[i]);
            if (table[i].dc <= 0) {
                cm_message_printf("Error: duty cycle %g is not positive, "
                                  "value replaced by 0.01.",
                                  table[i].dc);
                table[i].dc = 0.01;
            }
            else if (table[i].dc >= 1) {
                cm_message_printf("Error: duty cycle %g is 1 or larger, "
                                  "value replaced by 0.99.",
                                  table[i].dc);
                table[i].dc = 0.99;
            }
        }

        /* Allocate state data. */

        cm_event_alloc(0, sizeof (struct state));
        return;
    }

    table = STATIC_VAR(locdata);
    if (!table)
         return;
    state = (struct state *)cm_event_get_ptr(0, 0);

    if (CALL_TYPE != EVENT) {
        if (TIME == 0.0) {
            double phase;

            /* Set initial output and state data. */

            ctl = INPUT(cntl_in);
            ddc = get_dc(ctl, table, csize);

            phase = PARAM(init_phase);
            phase /= 360.0;
            if (phase < 0.0)
                phase += 1.0;

            /* When would a hypothetical previous transition have been? */

            state->last_time = delta * (1.0 - ddc - phase);
            if (state->last_time < 0.0) {
                state->last = ONE;
            } else {
                state->last = ZERO;
                state->last_time = -delta * phase;
            }
        }
        return;
    }

    /* Event call; either one requested previously or just before
     * a time-step is accepted.
     */

     if (TIME == 0.0) {
         OUTPUT_STATE(out) = state->last;
         OUTPUT_STRENGTH(out) = STRONG;
     }

     /* When is the next transition due? */

     ctl = INPUT(cntl_in);
     ddc = get_dc(ctl, table, csize);
     if (state->last)
         delta *= ddc;
     else
         delta *= (1.0 - ddc);
     when = state->last_time + delta;

     if (TIME >= when) {
         // If the frequency rose rapidly, the transition has been missed.
         // Force a shorter time-step and schedule then.

         if (!cm_analog_set_temp_bkpt(state->last_time + FACTOR2 * delta)) {
             OUTPUT_CHANGED(out) = FALSE;
             return;
         }

         /* Requested breakpoint was in the fixed past, or otherwise ignored:
          * request it later.
          */

         if (when > T(1) && !cm_analog_set_temp_bkpt((T(1) + when) / 2)) {
             OUTPUT_CHANGED(out) = FALSE;
             return;
         }

         /* Force output immediately. */

         when = TIME;
     }

     if (TIME >= state->last_time + FACTOR1 * delta) {
         /* TIME is reasonably close to transition time.  Request output. */

         state->last_time = when;
         state->last ^= ONE;
         OUTPUT_STATE(out) = state->last;
         OUTPUT_STRENGTH(out) = STRONG;
         OUTPUT_DELAY(out) = when - TIME;
         if (OUTPUT_DELAY(out) < 0.0)
             OUTPUT_DELAY(out) = 0.0;

         /* Request a call in the next half-cycle. */

	 if (state->last)
             delta *= ddc;
         else
             delta *= (1.0 - ddc);
         cm_event_queue(when + FACTOR2 * delta);
     } else {
         OUTPUT_CHANGED(out) = FALSE;

         if (TIME < state->last_time) {
             /* Output transition pending, nothing to do. */

             return;
         } else {
             /* Request a call nearer to transition time. */

             cm_event_queue(state->last_time + FACTOR2 * delta);
         }
     }
}
