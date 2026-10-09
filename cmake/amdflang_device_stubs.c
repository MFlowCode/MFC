/* AFAR 24.3.0's amdgcn libflang_rt.runtime.a references these (derived-type I/O
   tickets, verbose abort) but never defines them, so the device link fails. The
   stubs are weak and device-only: a real definition wins, and the host runtime
   is untouched. Calling one traps; MFC never does derived-type I/O on device. */
#pragma omp begin declare target device_type(nohost)

#define MFC_TRAP_STUB(name, sym)                                  \
    __attribute__((weak)) void name(void) __asm__(sym);           \
    __attribute__((weak)) void name(void) { __builtin_trap(); }

MFC_TRAP_STUB(mfc_stub_verbose_abort, "_Z22flang_rt_verbose_abortPKcz")
MFC_TRAP_STUB(mfc_stub_derived_in_continue,
              "_ZN7Fortran7runtime2io5descr15DerivedIoTicketILNS1_9DirectionE0EE8ContinueERNS0_9WorkQueueE")
MFC_TRAP_STUB(mfc_stub_derived_out_continue,
              "_ZN7Fortran7runtime2io5descr15DerivedIoTicketILNS1_9DirectionE1EE8ContinueERNS0_9WorkQueueE")
MFC_TRAP_STUB(mfc_stub_descriptor_in_begin,
              "_ZN7Fortran7runtime2io5descr18DescriptorIoTicketILNS1_9DirectionE0EE5BeginERNS0_9WorkQueueE")
MFC_TRAP_STUB(mfc_stub_descriptor_in_continue,
              "_ZN7Fortran7runtime2io5descr18DescriptorIoTicketILNS1_9DirectionE0EE8ContinueERNS0_9WorkQueueE")
MFC_TRAP_STUB(mfc_stub_descriptor_out_begin,
              "_ZN7Fortran7runtime2io5descr18DescriptorIoTicketILNS1_9DirectionE1EE5BeginERNS0_9WorkQueueE")
MFC_TRAP_STUB(mfc_stub_descriptor_out_continue,
              "_ZN7Fortran7runtime2io5descr18DescriptorIoTicketILNS1_9DirectionE1EE8ContinueERNS0_9WorkQueueE")

#pragma omp end declare target
