module event_handling
  use FoX_sax
contains

  subroutine characters_handler(chars)
    character(len=*), intent(in) :: chars

    print*, chars
  end subroutine
end module

program XMLreader
  use FoX_sax
  use event_handling
  type(xml_t) :: xp
  call open_xml_file(xp, 'input_3deq.xml')
  call parse(xp, characters_handler=characters_handler)
  call close_xml_t(xp)
end program
